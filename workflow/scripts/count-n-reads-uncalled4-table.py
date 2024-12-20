import polars as pl

# Define the path to your TSV file
file_path = "data/processed-data/R9/events-train.tsv"

unique_read_names = set()

reader = pl.read_csv_batched(
        file_path,
        separator="\t",
        batch_size=100000,
    )
batches = reader.next_batches(100)
counter = 0
while batches:
    # Update the set of unique read_name values
    df_current_batches = pl.concat(batches)
    read_names = df_current_batches["read_name"]
    unique_read_names.update(read_names)
    # DO Something
    batches = reader.next_batches(100)

unique_count = len(unique_read_names)

# Print the result
print(f"Number of unique read_name values: {unique_count}")
