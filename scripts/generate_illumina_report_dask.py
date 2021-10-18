import argparse
from illumina.report_generator_dask import prepare_report_using_dask

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--seqrun_id', required=True, help='Sequencing run id')
parser.add_argument('-j', '--stats_files_list', required=True, help='List of de-multiplexing Stats.json files')
parser.add_argument('-s', '--samplesheet_files', required=True, help='List of SampleSheet.csv files')
parser.add_argument('-o', '--output_file', required=True, help='output report path')
parser.add_argument('-t', '--template', required=True, help='Report template')
parser.add_argument('-w', '--workers', default=1, help='Number of workers for Dask LocalCluster')
parser.add_argument('-p', '--threads', default=1, help='Threads per Dask worker')
parser.add_argument('-m', '--memory_limit', default='1GB', help='Memory limit for Dask cluster')
parser.add_argument('-d', '--temp_dir', required=True, help='Dask temporary directory')
args = parser.parse_args()

seqrun_id = args.seqrun_id
stats_files_list = args.stats_files_list
samplesheet_files = args.samplesheet_files
output_file = args.output_file
template = args.template
workers = args.workers
threads = args.threads
memory_limit = args.memory_limit
temp_dir = args.temp_dir

if __name__=='__main__':
    prepare_report_using_dask(
        data_path=stats_files_list,
        samplesheets=samplesheet_files,
        seqrun_id=seqrun_id,
        template_path=template,
        n_workers=workers,
        threads_per_worker=threads,
        memory_limit=memory_limit,
        temp_dir=temp_dir,
        output_file=output_file)