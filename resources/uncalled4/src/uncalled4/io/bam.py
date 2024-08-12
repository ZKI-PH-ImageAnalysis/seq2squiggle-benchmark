import numpy as np
import pandas as pd
import sys
import pysam
import array
import os
import re
import time
import json
from collections import defaultdict, deque

from . import TrackIO

from ..moves import sam_to_ref_moves, sam_to_read_moves, INT32_NA
from ..tracks import AlnTrack, Alignment, AlnDF, LAYER_META, parse_layers
from ..ref_index import RefCoord
from .. import PoreModel
from .. import Config
from _uncalled4 import IntervalIndexI64, ValArrayI16, _RefCoord


REF_TAG  = "ur"
SAMP_TAG = "us"
LENS_TAG = "ul"
CURS_TAG = "uc"
STDS_TAG = "ud"

LAYER_PREFIXES = ["u","v","w","x","y","z"]
REQ_ALN_TAGS = [REF_TAG, LENS_TAG, CURS_TAG]

class BAM(TrackIO):
    FORMAT = "bam"

    def __init__(self, filename, write, tracks, track_count):
        TrackIO.__init__(self, filename, write, tracks, track_count)
        self.layer_tags = None

        if write:
            if not self.prms.buffered:
                if tracks.bam_in is None:
                    raise ValueError("No BAM template provided")
                self.header = tracks.bam_in.input.header.to_dict()
            else:
                self.header = self.prms.bam_header
            self.init_write_mode()
            self.input = None
        else:
            self.output = None
            self.init_read_mode()

    def init_read_mode(self):
        if self.conf.tracks.io.bam_in is None:
            raise ValueError("BAM output is only available if BAM input provided")

        self.track_in = None

        if self.prms.buffered:
            self.input = None
            self.header = self.prms.bam_header
        else:
            self.input = pysam.AlignmentFile(self.filename, "rb")

            self.header = self.input.header.to_dict()

        conf = None
        if "CO" in self.header:
            for line in self.header["CO"]:
                if not line.startswith("UNC:"): continue
                prms = json.loads(line[4:])
                
                if self.conf.tracks.ref_index is None:
                    self.conf.tracks.ref_index = prms["reference"]

                for name,vals in prms["tracks"].items():
                    c = self.conf.to_dict()
                    if not "pore_model" in c:
                        c["pore_model"] = {}
                    c["pore_model"].update(prms["models"][vals["model"]])
                    conf = Config(c)

                    if conf.read_index.paths is None:
                        conf.read_index.paths = vals["read"]["paths"]
                    else:
                        conf.read_index.paths = conf.read_index.paths + vals["read"]["paths"]

                    if conf.read_index.read_index is None or not os.path.exists(conf.read_index.read_index):
                        conf.read_index.read_index = vals["read"]["index"]

                    #TODO handle multiple tracks
                    self.track_in = self.init_track(self.input_name, vals["desc"], conf)
                    self.layer_tags = vals["layers"]

        if self.track_in is None:
            self.track_in = self.init_track(self.input_name, self.input_name, self.conf)

        self.load_moves = False
        for group,layer in parse_layers(self.conf.tracks.layers):
            if group == "moves":
                self.load_moves = True

        #if conf is None:
        #    conf = self.conf
        #self._init_tags(conf.tracks.io.bam_tags)
        #name = os.path.basename(self.filename)
        #self.track_in = self.init_track(name, name, conf)

        self.read_id_index = None

        self.in_alns = None

    def reset(self):
        self.input.reset()

    def iter_read(self, read_id):
        if self.read_id_index is None:
            self.read_id_index = pysam.IndexedReads(self.input)
            self.read_id_index.build()
        ret = self.read_id_index.find(read_id)
        return ret

    #def _init_tags(self, tags):
    #    self.tags = dict()
    #    for t in tags:
    #        label, tag = t.split(":")
    #        self.tags[label] = tag

    def init_write_mode(self):
        if self.conf.tracks.io.bam_out is None:
            self.output = None
            return

        #self._init_tags(self.prms.bam_tags)

        TrackIO.init_write_mode(self)

        #TODO load from AlnTrack instance (initialized by Tracks)
        self.model = PoreModel(params=self.conf.pore_model)
        #self.kmer_trim = self.model.kmer_trim

        if self.prms.buffered:
            self.out_buffer = list()
            self.output = None
            return

        #Store config toml in single-line comment with newlines replaced by semicolons
        conf_line = self.conf.to_toml() \
                        .replace(";", "\\;") \
                        .replace("\n", ";")   

        params = {
            "tracks" : {
                self.track_out.name : {
                    "desc" : self.track_out.desc,
                    "model" : self.track_out.model.name, #self.model.name, #
                    "read" : {"paths" : self.conf.read_index.paths, "index" : self.conf.read_index.read_index},
                    "layers" : {                                                                 
                        LENS_TAG : {"name" : "dtw.length"},                                      
                        CURS_TAG : {"name" : "dtw.current", "scale" : self.track_out.model.inorm_scale},   
                        STDS_TAG : {"name" : "dtw.current_sd", "scale" : self.track_out.model.inorm_scale},
            }}},
            "models" : {
                self.track_out.model.name : self.track_out.model.params_to_dict()
            },
            "reference" : self.conf.tracks.ref_index
        }

        if not "CO" in self.header:
            self.header["CO"] = list()
        #self.header["CO"].append("UNC:" + conf_line)
        self.header["CO"].append("UNC:" + json.dumps(params))

        self.output = pysam.AlignmentFile(self.conf.tracks.io.bam_out, "wb", header=self.header)#template=self.input)

    def get_alns(self, read_id):
        self._init_alns()
        return self.in_alns[read_id]

    def get_aln(self, read_id, ref_name, ref_start):
        self._init_alns()
        for aln in self.in_alns[read_id]:
            if aln.reference_name == ref_name and aln.reference_start == ref_start:
                return aln
        return None

    MIN_I16 = np.iinfo(np.int16).min+1
    MAX_I16 = np.iinfo(np.int16).max
    NA_I16 = ValArrayI16.NA

    def _init_alns(self):
        if self.in_alns is None:
            self.in_alns = defaultdict(list)
            for aln in self.iter_sam():
                self.in_alns[aln.query_name].append(aln)
            self.input.reset()

    LAYER_TAGS = {
        "dtw.length"     : (LENS_TAG, "h"),
        "dtw.current"    : (CURS_TAG, "h"),
        "dtw.current_sd" : (STDS_TAG, "h"),
        "dtwcmp.dist"    : (None, "f"),
    }

    def _current_to_tag(self, tag, vals):
        if len(vals) == 0:
            return False

        ints = np.round(np.array(vals) * self.model.current.INORM_SCALE)
        ints[np.isnan(vals)] = self.NA_I16
        arr = array.array("h", ints.astype(np.int16))
        self.bam.set_tag(tag, arr)

        return True

    def write_alignment(self, aln):
        self.bam = aln.sam

        refs = aln.seq.coord.bounds

        self._current_to_tag(CURS_TAG, aln.dtw.current)
        self._current_to_tag(STDS_TAG, aln.dtw.current_sd)

        start_pad = list()
        start = -aln.dtw.samples.start
        while start < self.MIN_I16:
            start_pad.append(self.MIN_I16)
            start -= self.MIN_I16
        if start < 0:
            start_pad.append(start)

        lens = np.concatenate([start_pad, aln.dtw.samples.to_runlen()])#.astype(np.int16)
        #old = lens
        lens = lens.astype(np.int16)

        #mis = old != lens

        self.bam.set_tag(REF_TAG, array.array("i", refs))
        self.bam.set_tag(LENS_TAG, array.array("h", lens))

        if self.prms.bam_f5c:
            rst,ren = refs[0], refs[-1]-self.model.k+1
            if self.model.reverse:
                rst,ren = ren,rst
            si_tag = (aln.dtw.samples.start, aln.dtw.samples.end, rst,ren)
            ss_tag = list()
            i = 0
            prev = None
            ref_gap = 0
            for j in range(0,len(refs), 2):
                st,en = refs[j:j+2]
                if j+3 >= len(refs):
                    en -= self.model.K - 1
                for r in range(st,en):
                    intv = aln.dtw.samples.get_interval(i)
                    c = aln.dtw.current[i]
                    if not intv.is_valid() or (prev is not None and prev.start == intv.start):
                        ref_gap += 1
                    elif not np.isnan(c) and intv.is_valid():
                        if ref_gap > 0:
                            ss_tag.append(f"{ref_gap}D")
                            ref_gap = 0

                        if prev is not None and intv.start > prev.end:
                            l = intv.start - prev.end
                            ss_tag.append(f"{l}I")

                        ss_tag.append(f"{len(intv)},")
                        prev = intv
                        
                    else:
                        sys.stderr.write(f"Error writing ss tag for {aln.read.id} ({int(intv.is_valid())}{int(np.isnan(c))})\n")
                    i += 1
                    
                if j+3 < len(refs):
                    ref_gap += refs[j+2] - refs[j+1]
                    #ss_tag.append(f"{l}D")

            sc,sh = aln.get_scaled_norm(self.model.pa_mean, self.model.pa_stdv)
            self.bam.set_tag("sc", 1/sc)
            self.bam.set_tag("sh", -sh*sc)

            self.bam.set_tag("ss", "".join(ss_tag))
            self.bam.set_tag("si", ",".join(map(str,si_tag)))

        if self.prms.buffered:
            self.out_buffer.append(self.bam.to_dict())
        else:
            self.output.write(self.bam)

    def write_buffer(self, bams):
        header = pysam.AlignmentHeader.from_dict(self.header)
        for b in bams:
            bam = pysam.AlignedSegment.from_dict(b, header)
            #bam = pysam.AlignedSegment.fromstring(b, header)
            self.output.write(bam)

    def iter_str_chunks(self, group_seqs=False):
        read_ids = set()
        bams = list()
        prev_seq = None
        for bam in self.iter_sam():
            read_ids.add(bam.query_name)
            if bam.has_tag("pi"): # splitted re
                read_ids.add(bam.get_tag("pi"))
            bams.append(bam.to_string())

            if len(bams) >= self.prms.bam_chunksize and (not group_seqs or bam.reference_name != prev_seq):
                yield(read_ids, bams)
                read_ids = set()
                bams = list()

            prev_seq = bam.reference_name

        if len(bams) > 0:
            yield(read_ids, bams)

    #TODO more query options
    def iter_sam(self, unmapped=False):
        #if self.conf.tracks.ref_bounds is not None
        if self.conf.tracks.ref_bounds is None:
            itr = self.input
        else:
            b = self.conf.tracks.ref_bounds
            if b.has_bounds:
                itr = self.input.fetch(b.name, b.start, b.end)
            else:
                itr = self.input.fetch(b.name)

        if unmapped:
            mapped = lambda a: True
        else:
            mapped = lambda a: not a.is_unmapped

        read_filter = self.tracks.read_index.read_filter
        if read_filter is not None:
            filt = lambda a: a.query_name in read_filter
        else:
            filt = lambda a: True

        valid = lambda a: mapped(a) and filt(a)
            
        if self.conf.tracks.max_reads is None:
            for a in itr:
                if valid(a): yield a
        else:
            n = 0
            for a in itr:
                if valid(a): 
                    yield a
                    n += 1
                if n == self.conf.tracks.max_reads:
                    break

    def iter_alns(self):#, layers=None, track_id=None, coords=None, aln_id=None, read_id=None, fwd=None, full_overlap=None, ref_index=None):
        for sam in self.iter_sam():
            aln = self.sam_to_aln(sam)
            yield aln

    def _tag_to_layer(self, sam, tag):
        if sam.has_tag(tag):
            return None
        if not tag in self.layer_tags:
            raise ValueError(f"Unknown layer tag: {tag}")
        l = self.layer_tags[tag]
        vals = np.array(sam.get_tag(tag))
        na = vals == self.NA_I16
        vals = (vals - l["shift"]) / l["scale"]
        vals[na] = pd.nan
        return (l["name"], vals)

    def sam_to_aln(self, sam, load_moves=None):
        if load_moves is None: load_moves = self.load_moves

        sys.stdout.flush()
        has_dtw = True
        for tag in REQ_ALN_TAGS:
            if not sam.has_tag(tag):
                has_dtw = False
                break

        aln = None

        if sam.has_tag("pi"): # splitted read, should have "sp".
            read = self.tracks.read_index[sam.get_tag("pi")]
            read = self.tracks.read_index.process_splitted_read(read, sam)
        else:
            read = self.tracks.read_index[sam.query_name] #None)

        if not has_dtw and read is None:
            return None

        fwd = int(not sam.is_reverse)

        if has_dtw:
            layers = dict()

            for tag, meta in self.layer_tags.items():
                if not sam.has_tag(tag):
                    continue
                vals = np.array(sam.get_tag(tag))
                na = vals == self.NA_I16
                vals = (vals - meta.get("shift",0)) / meta.get("scale", 1)
                vals[na] = np.nan
                layers[meta["name"]] = vals

            refs = sam.get_tag(REF_TAG)

            #coords = RefCoord(sam.reference_name, refs, fwd)
            #aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read, sam.reference_id, coords, sam)
            aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read, sam, refs)

            length = layers["dtw.length"]
            mask = length >= 0
            layers["dtw.start"] = np.pad(np.cumsum(np.abs(length)), (1,0))[:-1][mask]

            length = length[mask]

            try:
                lna = (length == 0) & np.isnan(layers["dtw.current"])
            except ValueError:
                return None

            length[lna] = -1
            layers["dtw.length"] = length.astype(np.int32)

            dtw = AlnDF(aln.seq, layers["dtw.start"], layers["dtw.length"], layers.get("dtw.current", None), layers.get("dtw.current_sd", None))#start, length, current, stdv)

            aln.set_dtw(dtw)

        moves = None
        if load_moves:

            if self.tracks.index is not None:
                moves = sam_to_ref_moves(self.conf, self.tracks.index, read, sam)
            else:
                moves = sam_to_read_moves(read, sam)
                #moves = moves.slice(1,len(moves))

            if moves is not None:
                if aln is None:
                    #coords = RefCoord(sam.reference_name, moves.index, fwd)
                    #aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read, sam.reference_id, coords, sam)
                    aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read, sam, moves.index)
                else:
                    i = max(0, moves.index.start - aln.seq.index.start)
                    j = min(len(moves), len(moves) + moves.index.end -  aln.seq.index.end)
                    moves = moves.slice(i,j)
                    
                aln.set_moves(moves)

        if moves is not None and has_dtw:
            aln.calc_mvcmp()
        elif aln is None:
            #coords = RefCoord(sam.reference_name, sam.reference_start, sam.reference_end, fwd)
            #aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read, sam.reference_id, coords, sam)
            aln = self.tracks.init_alignment(self.track_in.name, self.next_aln_id(), read, sam, sam.reference_start, sam.reference_end)
        
        return aln

    def _parse_sam(self, sam, layer_names):
        aln = self.sam_to_aln(sam)

        layers = aln.to_pandas(layer_names, index=["aln.track", "seq.fwd", "seq.pos", "aln.id"])
        #layers["track.name"] = self.track_in.name
        #layers = layers.set_index("track.name", append=True)\
        #               .reorder_levels(["track.name", "seq.fwd", "seq.pos", "aln.id"])

        attr = aln.attrs()
        aln_df = pd.DataFrame([attr], columns=attr._fields)
        #aln_df["track.name"] = self.track_in.name
        aln_df.set_index(["aln.track","id"], inplace=True)

        return aln_df, layers

    def query(self, layers, coords, index=["aln.id","seq.pos"], read_id=None, full_overlap=False):
        itr = self.input.fetch(coords.name, coords.start, coords.end)

        track_alns = defaultdict(list)
        track_layers = defaultdict(list)

        if read_id is not None and len(read_id) > 0:
            read_ids = set(read_id)
        else:
            read_ids = None

        for sam in itr:
            if read_ids is not None and sam.query_name is not read_ids:
                continue
            aln = self.sam_to_aln(sam)

            track_alns[self.track_in.name].append(aln.attrs())

            l = aln.to_pandas(layers, index=index, bounds=coords)
            track_layers[self.track_in.name].append(l)

        track_alns = {
            t : pd.DataFrame(l, columns=Alignment.Attrs._fields)\
                  .set_index("id")
            for t,l in track_alns.items()}
        track_layers = {t : pd.concat(l) for t,l in track_layers.items()}
            
        return track_alns, track_layers
    
    def iter_bed(self, fname):
        coords = pd.read_csv(fname, sep="\t", names=["name","start","end"])
        coords["pac_start"] = [self.tracks.index.get_pac_offset(r) for r in coords["name"]]
        coords = coords.sort_values("pac_start")
        #for chunk in pd.read_csv(fname, sep="\t", names=["name","start","end"], chunksize=100):
        for i,b in coords.iterrows():
            for aln in self.input.fetch(b["name"], b["start"], b["end"]):
                yield aln

    def iter_refs(self, layers, coords=None, aln_id=None, read_id=None, fwd=None, chunksize=None, full_overlap=False):

        if coords is not None:
            itr = self.input.fetch(coords.ref_name, coords.refs.min(), coords.refs.max())
        elif self.prms.bed_filter is not None:
            itr = self.iter_bed(self.prms.bed_filter)
        else:
            itr = self.input.fetch()

        if fwd is not None:
            strands = [int(fwd)]
        else:
            strands = [1, 0]
        
        new_layers = list()
        new_alns = list()

        layer_df = pd.DataFrame()
        aln_df = pd.DataFrame()

        prev_ref = None
        prev_start = 0
        #########prev_aln = 0

        t = time.time()
        for sam in itr:
            aln = self.sam_to_aln(sam)

            next_aln = aln.attrs()
            next_layers = aln.to_pandas(layers, index=["seq.pac", "seq.fwd", "aln.id"])
            next_ref = next_aln.ref_name

            next_start = sam.reference_start #next_layers.index.get_level_values("seq.pac").min()

            first = prev_ref is not None
            new = next_ref != prev_ref
            loong = next_start - prev_start > 0
            big = len(new_alns) > self.prms.aln_chunksize
            if first and (new or (loong and big)):
            #if prev_ref is not None and (next_ref != prev_ref or (next_start - prev_start > 0 and len(new_alns) > self.prms.aln_chunksize)):

                aln_df = pd.concat([aln_df, pd.DataFrame(new_alns, columns=new_alns[0]._fields).set_index("id")])
                layer_df = pd.concat([layer_df] + new_layers).sort_index()
                new_alns = list()
                new_layers = list()
                
                #ret_layers = layer_df.loc[(slice(None),slice(prev_start,next_start-1)), :]
                #ret_alns = aln_df.loc[ret_layers.index.droplevel(["seq.fwd", "seq.pac"]).unique()]
                if next_ref == prev_ref:
                    ret_layers = layer_df.loc[prev_start:next_start-1]
                    ret_alns = aln_df.loc[ret_layers.index.unique("aln.id")]

                    layer_df = layer_df.loc[next_start:]
                    aln_df = aln_df.loc[layer_df.index.unique("aln.id")]
                else:
                    ret_layers = layer_df#.loc[prev_start:]
                    ret_alns = aln_df#.loc[ret_layers.index.unique("aln.id")]
                    layer_df = layer_df.iloc[:0]
                    aln_df = aln_df.iloc[:0]

                t = time.time()

                if len(ret_alns) > 0:
                    yield (ret_alns, ret_layers)

                prev_ref = next_ref
                prev_start = next_start

            elif prev_ref != next_ref:
                prev_ref = next_ref
                prev_start = next_start

            new_alns.append(next_aln)
            new_layers.append(next_layers)

        new_alns = pd.DataFrame(new_alns, columns=new_alns[0]._fields).set_index("id")
        aln_df = pd.concat([aln_df] + [new_alns]).sort_index()
        layer_df = pd.concat([layer_df] + new_layers).sort_index()

        ret_layers = layer_df.loc[prev_start:]
        ret_alns = aln_df.loc[ret_layers.index.droplevel(["seq.fwd", "seq.pac"]).unique()]

        if len(ret_alns) > 0:
            yield (ret_alns, ret_layers)



    def close(self):
        if self.input is not None:
            self.input.close()
        if self.output is not None:
            self.output.close()

