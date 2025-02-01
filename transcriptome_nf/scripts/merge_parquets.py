
import os, sys
import pyarrow as pa
import pyarrow.parquet as pq

def merge_parquet_files(files, output_file):
    #files = [os.path.join(path,f) for f in os.listdir(path)]
    schema = pq.ParquetFile(files[0]).schema_arrow
    with pq.ParquetWriter(output_file, schema=schema) as writer:
        for file in files:
            writer.write_table(pq.read_table(file, schema=schema))

if __name__ == "__main__":
    # e.g. python merge_parquets.py file1.parquet file2.parquet ... merged.parquet
    #parquet_path, merged_file = sys.argv[1:]
    #merge_parquet_files(parquet_path, merged_file)
    *parquet_files, merged_file = sys.argv[1:]
    merge_parquet_files(parquet_files, merged_file)
