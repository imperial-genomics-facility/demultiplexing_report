import argparse, logging
from illumina.report_generator_for_db import create_plot_json_for_database

parser = argparse.ArgumentParser()
parser.add_argument('-i','--seqrun_id', required=True, help='Sequencing run id')
parser.add_argument('-j','--stats_files_list', required=True, help='List of de-multiplexing Stats.json files')
parser.add_argument('-s','--samplesheet_files', required=True, help='List of SampleSheet.csv files')
parser.add_argument('-o','--output_dir', required=True, help='output report dir')
parser.add_argument('-t','--samplesheet_tag', required=True, help='Samplesheet tag')
args = parser.parse_args()

seqrun_id = args.seqrun_id
stats_files_list = args.stats_files_list
samplesheet_files = args.samplesheet_files
output_dir = args.output_dir
samplesheet_tag = args.samplesheet_tag

if __name__=='__main__':
    try:
        create_plot_json_for_database(
            run_name=seqrun_id,
            samplesheet_tag=samplesheet_tag,
            stat_files=stats_files_list,
            samplesheet_files=samplesheet_files,
            output_dir=output_dir)
    except Exception as e:
        logging.error(e)