import unittest, json
import pandas as pd
import sys
#sys.path.append('/home/vmuser/github/demultiplexing_report')
from illumina.report_generator import read_bcl2fastq_stats_data_from_pandas
from illumina.report_generator import read_data_via_pandas
from illumina.report_generator import get_stats_summary_table
from illumina.report_generator import get_samplesheet_records
from illumina.report_generator import get_flowcell_summary_plots
from illumina.report_generator import get_per_lane_sample_dist_plot
from illumina.report_generator import get_project_summary_html_table
from illumina.report_generator import get_demult_per_lane_demult_table_data
from illumina.report_generator import get_flowcell_project_summary_plot
from illumina.report_generator import get_undetermined_plot
from illumina.report_generator import get_undetermined_table
from illumina.report_generator import create_demux_html_report
from illumina.report_generator import combine_data_and_create_report
from illumina.report_generator import prepare_report_using_pandas


class Report_generator1(unittest.TestCase):
    def setUp(self):
        #self.stats_json = '/home/vmuser/github/demultiplexing_report/data/illumina/Stats_v3_1.json'
        #self.samplesheet = '/home/vmuser/github/demultiplexing_report/data/illumina/SampleSheet_v3_1.csv'
        self.stats_json = 'data/illumina/Stats_v3_1.json'
        self.samplesheet = 'data/illumina/SampleSheet_v3_1.csv'

    def tearDown(self):
        pass

    def test_read_bcl2fastq_stats_data_from_pandas(self):
        with open(self.stats_json, 'r') as jp:
            data = json.load(jp)
            row_l, row_s, unknown_df = \
                read_bcl2fastq_stats_data_from_pandas(data)
            self.assertEqual(len(row_l), 1)
            self.assertEqual(len(row_s), 95)
            self.assertEqual(len(unknown_df), 1000)
            row_l = pd.DataFrame(row_l)
            self.assertEqual(row_l['Total_cluster_raw'].values.tolist()[0], 8624790)
            self.assertEqual(row_l['Total_cluster_pf'].values.tolist()[0], 3150921)
            row_s = pd.DataFrame(row_s)
            self.assertEqual(
                row_s[row_s['Sample_Name']=='INS0008D8']['Num_reads'].values[0], 29157)
            unknown_df = pd.DataFrame(unknown_df)
            self.assertEqual(
                unknown_df[unknown_df['Barcode']=='GGCAAGTT+CCCCCCCC']['Reads'].values[0], 440)

    def test_read_data_via_pandas(self):
        summary_records, sample_records, undetermined_records = \
            read_data_via_pandas(data_path=[self.stats_json])
        self.assertEqual(
            summary_records['Total_cluster_raw'].values.tolist()[0], 8624790)
        self.assertEqual(
            summary_records['Total_cluster_pf'].values.tolist()[0], 3150921)
        self.assertEqual(
                sample_records[sample_records['Sample_Name']=='INS0008D8']['Num_reads'].values[0], 29157)
        self.assertEqual(
                undetermined_records[undetermined_records['Barcode']=='GGCAAGTT+CCCCCCCC']['Reads'].values[0], 440)

    def test_get_stats_summary_table(self):
        summary_records, sample_records, undetermined_records = \
            read_data_via_pandas(data_path=[self.stats_json])
        records = \
            get_stats_summary_table(
                sum_df=summary_records,
                lane_sample_df=sample_records)
        self.assertEqual(
            records[records['Sample_Name']=='INS0008D8']['PF Clusters'].values[0], '29,157')

    def test_get_samplesheet_records(self):
        records = \
            get_samplesheet_records(
                samplesheets=[self.samplesheet])
        self.assertEqual(
            records[records['Sample_ID']=='IGF118539']['Sample_Project'].values[0], 'IGFQ001220')
        self.assertEqual(
            len(records.index), 95)

    def test_get_flowcell_summary_plots(self):
        summary_records, sample_records, undetermined_records = \
            read_data_via_pandas(data_path=[self.stats_json])
        summary_plot1, summary_plot2 = \
            get_flowcell_summary_plots(
                summary_data=summary_records)
        self.assertTrue(isinstance(summary_plot1, str))
        self.assertTrue(isinstance(summary_plot2, str))

    def get_per_lane_sample_dist_plot(self):
        bg_colors = ['rgba(75, 192, 192, 0.2)']
        border_colors = ['rgba(75, 192, 192, 1)']
        _, sample_records, _ = \
            read_data_via_pandas(data_path=[self.stats_json])
        lane_plots, plot_height = \
            get_per_lane_sample_dist_plot(
                sample_data=sample_records,
                bg_colors=bg_colors,
                border_colors=border_colors)
        self.assertTrue(isinstance(lane_plots, dict))
        self.assertEqual(plot_height, 1000)

if __name__ == '__main__':
  unittest.main()