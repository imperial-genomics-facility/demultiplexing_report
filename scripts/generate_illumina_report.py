import argparse
from ..illumina.report_generator import prepare_report_using_pandas

parser = argparse.ArgumentParser()
parser.add_argument('-i','--seqrun_id', required=True, help='Sequencing run id')
parser.add_argument('-j','--stats_files_list', required=True, help='List of de-multiplexing Stats.json files')
parser.add_argument('-s','--samplesheet_files', required=True, help='List of SampleSheet.csv files')
parser.add_argument('-o','--output_file', required=True, help='output report path')
parser.add_argument('-t','--template', required=True, help='Report template')
args = parser.parse_args()

seqrun_id = args.seqrun_id
stats_files_list = args.stats_files_list
samplesheet_files = args.samplesheet_files
output_file = args.output_file
template = args.template

if __name__=='__main__':
    prepare_report_using_pandas(
        data_path=stats_files_list,
        samplesheets=samplesheet_files,
        seqrun_id=seqrun_id,
        template_path=template,
        output_file=output_file)