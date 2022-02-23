#!/usr/bin/env python
desc= "Concatenate file(s) containing median current intensity values per position values."

# Import required libraries:
import argparse
import os

import duckdb

def process_contig_group(files, output_file, num_threads, countig_groups_fn, group_id):
    
    def parse_contig_groups(tsv_fn, group_id):
        my_group_contigs = []
        with open(tsv_fn) as f:
            for line in f:
                line = line.strip()
                contig, group_name = line.split()
                if group_name == group_id:
                    my_group_contigs.append(contig)
        
        assert len(my_group_contigs) > 0
        return my_group_contigs

    # set up DuckDB    
    duckdb_pragma_1 = "PRAGMA enable_object_cache"
    duckdb_pragma_2 = f"PRAGMA threads={num_threads}"

    con = duckdb.connect(database=':memory:')
    con.execute(duckdb_pragma_1)
    con.execute(duckdb_pragma_2)
    
    group_parquet_glob = f"*_dset/contig_group={group_id}/*parquet"
    print(group_parquet_glob)
    sql_option_1 = f"CREATE VIEW contig_grp AS SELECT * FROM parquet_scan('{group_parquet_glob}')"
    
    sql_option_2 = f"CREATE TABLE contig_grp AS SELECT * FROM parquet_scan('{group_parquet_glob}')" 
    if debug_me == True:
        print(option_1)
        print(option_1)
    con.execute(f"{sql_option_2}")

    # optional, needs to be benchmarked
    sql_index = "CREATE INDEX con_pos_kmer ON contig_grp (contig, position, reference_kmer)"
    con.execute(f"{sql_index}")

    contig_list =  parse_contig_groups(countig_groups_fn, group_id)

    for contig_id in contig_list:
        out_csv_fn = f"{group_id}_{contig_id}_merged.csv"
        contig_select_sql = f"COPY (SELECT contig, position, reference_kmer, median(median), sum(coverage)  \
            FROM contig_grp \
            WHERE contig='{contig_id}' \
            GROUP BY contig, position, reference_kmer) TO '{out_csv_fn}' WITH (HEADER 1, DELIMITER ',')"
        print(contig_select_sql)
        con.execute(f"{contig_select_sql}")
        command_pigz = f"pigz -p{num_threads} {out_csv_fn}"
       
        os.system(f"{command_pigz}")
    



def main():
    parser  = argparse.ArgumentParser(description=desc)

    parser.add_argument('-d', '--dir', nargs='+' ,help='Input top dir to process.')
    parser.add_argument('-o', '--output', help='Output dir')
    parser.add_argument('-t', '--threads', help='threads used by DuckDB')
    parser.add_argument('-l', '--countig_groups', help='TSV file contig_name contig_group')
    parser.add_argument('-g', '--group_id', help='contig_group id to process')

    a = parser.parse_args()
    
    #Read, parse and merge data from individual eventalign files:
    process_contig_group(a.input, a.output, a.threads, a.countig_groups, a.group_id)
    
if __name__=='__main__':
    debug_me = False
    main()
