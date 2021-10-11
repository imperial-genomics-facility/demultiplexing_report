import pandas as pd
import unittest, json, sys, os
#sys.path.append('/home/vmuser/github/demultiplexing_report')
from illumina.report_generator import read_bcl2fastq_stats_data_from_pandas
from illumina.report_generator import read_data_via_pandas
from illumina.report_generator import get_stats_summary_table
from illumina.report_generator import get_samplesheet_records
from illumina.report_generator import get_flowcell_summary_plots
from illumina.report_generator import get_per_lane_sample_dist_plot
from illumina.report_generator import get_project_summary_html_table
from illumina.report_generator import get_per_lane_demult_table_data
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
        summary_records, _, _ = \
            read_data_via_pandas(data_path=[self.stats_json])
        summary_plot1, summary_plot2 = \
            get_flowcell_summary_plots(
                summary_data=summary_records)
        self.assertTrue(isinstance(summary_plot1, str))
        self.assertTrue(isinstance(summary_plot2, str))

    def test_get_per_lane_sample_dist_plot(self):
        bg_colors = ['rgba(75, 192, 192, 0.2)']
        border_colors = ['rgba(75, 192, 192, 1)']
        sum_df, sample_records, _ = \
            read_data_via_pandas(data_path=[self.stats_json])
        all_samplesheet_data = \
            get_samplesheet_records([self.samplesheet])
        sample_data = \
            get_stats_summary_table(
                sum_df,
                sample_records)
        merged_sample_data = \
            sample_data.set_index('Sample_ID').\
            join(
                all_samplesheet_data.set_index('Sample_ID')['Sample_Project'],
                how='left').\
            reset_index().\
            drop_duplicates()
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].\
                map(lambda x: x.replace(',',''))
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].astype(int)
        lane_plots, plot_height = \
            get_per_lane_sample_dist_plot(
                sample_data=merged_sample_data,
                bg_colors=bg_colors,
                border_colors=border_colors)
        self.assertTrue(isinstance(lane_plots, dict))
        self.assertEqual(plot_height, 1000)

    def test_get_project_summary_html_table(self):
        sum_df, sample_records, _ = \
            read_data_via_pandas(data_path=[self.stats_json])
        all_samplesheet_data = \
            get_samplesheet_records([self.samplesheet])
        sample_data = \
            get_stats_summary_table(
                sum_df,
                sample_records)
        merged_sample_data = \
            sample_data.set_index('Sample_ID').\
            join(
                all_samplesheet_data.set_index('Sample_ID')['Sample_Project'],
                how='left').\
            reset_index().\
            drop_duplicates()
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].\
                map(lambda x: x.replace(',',''))
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].astype(int)
        project_summary_data = \
            get_project_summary_html_table(
                sample_data=merged_sample_data)
        self.assertEqual(type(project_summary_data), str)
        self.assertTrue('<td>IGFQ001220</td>' in project_summary_data)

    def test_get_per_lane_demult_table_data(self):
        sum_df, sample_records, _ = \
            read_data_via_pandas(data_path=[self.stats_json])
        all_samplesheet_data = \
            get_samplesheet_records([self.samplesheet])
        sample_data = \
            get_stats_summary_table(
                sum_df,
                sample_records)
        merged_sample_data = \
            sample_data.set_index('Sample_ID').\
            join(
                all_samplesheet_data.set_index('Sample_ID')['Sample_Project'],
                how='left').\
            reset_index().\
            drop_duplicates()
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].\
                map(lambda x: x.replace(',',''))
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].astype(int)
        table_data = \
            get_per_lane_demult_table_data(
                sample_data=merged_sample_data)
        self.assertTrue(isinstance(table_data, dict))
        self.assertEqual(len(table_data.keys()), 1)
        self.assertTrue(1 in table_data.keys())

    def test_get_flowcell_project_summary_plot(self):
        sum_df, sample_records, _ = \
            read_data_via_pandas(data_path=[self.stats_json])
        all_samplesheet_data = \
            get_samplesheet_records([self.samplesheet])
        sample_data = \
            get_stats_summary_table(
                sum_df,
                sample_records)
        merged_sample_data = \
            sample_data.set_index('Sample_ID').\
            join(
                all_samplesheet_data.set_index('Sample_ID')['Sample_Project'],
                how='left').\
            reset_index().\
            drop_duplicates()
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].\
                map(lambda x: x.replace(',',''))
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].astype(int)
        project_summary_plot = \
            get_flowcell_project_summary_plot(
                summary_data=sum_df,
                sample_data=merged_sample_data)
        self.assertTrue(isinstance(project_summary_plot, str))

    def test_get_undetermined_plot(self):
        _, _, undetermined_data = \
            read_data_via_pandas(data_path=[self.stats_json])
        bg_colors = ['rgba(75, 192, 192, 0.2)']
        border_colors = ['rgba(75, 192, 192, 1)']
        undetermined_plots = \
            get_undetermined_plot(
                undetermined_data=undetermined_data,
                bg_colors=bg_colors,
                border_colors=border_colors)
        self.assertTrue(isinstance(undetermined_plots, dict))
        self.assertTrue(1 in undetermined_plots.keys())

    def test_get_undetermined_table(self):
        _, _, undetermined_data = \
            read_data_via_pandas(data_path=[self.stats_json])
        undetermined_tables = \
            get_undetermined_table(
                undetermined_data=undetermined_data)
        self.assertTrue(isinstance(undetermined_tables, dict))
        self.assertTrue(1 in undetermined_tables.keys())

    def test_create_demux_html_report(self):
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_pandas(data_path=[self.stats_json])
        if os.path.exists('/tmp/t1.html'):
            os.remove('/tmp/t1.html')
        combine_data_and_create_report(
            sum_df=sum_df,
            lane_sample_df=lane_sample_df,
            undetermined_data=undetermined_data,
            samplesheets=[self.samplesheet],
            seqrun_id='TEST 1',
            template_path='template/illumina_report_v1.html',
            output_file='/tmp/t1.html')
        self.assertTrue(os.path.exists('/tmp/t1.html'))

    def test_prepare_report_using_pandas(self):
        samplesheet_list = '/tmp/samplesheet_list.csv'
        stats_json_list = '/tmp/stats_json.csv'
        if os.path.exists('/tmp/t1.html'):
            os.remove('/tmp/t1.html')
        with open(samplesheet_list, 'w') as fp:
            fp.write(self.samplesheet)
        with open(stats_json_list, 'w') as fp:
            fp.write(self.stats_json)
        prepare_report_using_pandas(
            data_path=stats_json_list,
            samplesheets=samplesheet_list,
            seqrun_id='TEST 1',
            template_path='template/illumina_report_v1.html',
            output_file='/tmp/t1.html')
        self.assertTrue(os.path.exists('/tmp/t1.html'))

if __name__ == '__main__':
  unittest.main()