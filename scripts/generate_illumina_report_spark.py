import argparse
from illumina.report_generator_spark import prepare_report_using_spark

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--seqrun_id', required=True, help='Sequencing run id')
parser.add_argument('-j', '--stats_files_list', required=True, help='List of de-multiplexing Stats.json files')
parser.add_argument('-s', '--samplesheet_files', required=True, help='List of SampleSheet.csv files')
parser.add_argument('-o', '--output_file', required=True, help='output report path')
parser.add_argument('-t', '--template', required=True, help='Report template')
parser.add_argument('-p', '--threads', default=1, help='Number of threads for Spark cluster')
parser.add_argument('-m', '--memory_limit', default=4, help='Memory limit for Spark executor')
args = parser.parse_args()

seqrun_id = args.seqrun_id
stats_files_list = args.stats_files_list
samplesheet_files = args.samplesheet_files
output_file = args.output_file
template = args.template
threads = args.threads
memory_limit = args.memory_limit

if __name__=='__main__':
    prepare_report_using_spark(
        data_path=stats_files_list,
        samplesheets=samplesheet_files,
        seqrun_id=seqrun_id,
        template_path=template,
        threads=threads,
        executor_memory_in_gb=memory_limit,
        output_file=output_file)