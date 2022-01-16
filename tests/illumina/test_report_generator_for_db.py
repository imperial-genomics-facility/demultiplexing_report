import unittest
from illumina.report_generator_for_db import create_plot_json_for_database
from illumina.report_generator_for_db import get_flowcell_summary_data
from illumina.report_generator_for_db import get_per_lane_sample_dist_plot
from illumina.report_generator_for_db import get_project_summary_html_table
from illumina.report_generator_for_db import get_per_lane_demult_table_data
from illumina.report_generator_for_db import get_flowcell_project_summary_plot_for_db
from illumina.report_generator_for_db import get_undetermined_table
from illumina.report_generator_for_db import get_undetermined_plot
from illumina.report_generator import (
    read_data_via_pandas,
    get_stats_summary_table,
    get_samplesheet_records)

class Test_report_generator_for_db(unittest.TestCase):
    def setUp(self):
        self.stats_json = 'data/illumina/Stats_v3_1.json'
        self.samplesheet = 'data/illumina/SampleSheet_v3_1.csv'
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_pandas(
                data_path=[self.stats_json])
        sample_data = \
            get_stats_summary_table(
                sum_df,
                lane_sample_df)
        all_samplesheet_data = \
            get_samplesheet_records([self.samplesheet])
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
        self.sum_df = sum_df
        self.merged_sample_data = merged_sample_data
        self.undetermined_data = undetermined_data
    def tearDown(self):
        pass

    def test_get_flowcell_summary_data(self):
        (labels, total_cluster_raw, total_cluster_pf, total_yield) = \
            get_flowcell_summary_data(summary_data=self.sum_df)
        self.assertTrue('Lane 1' in labels)
        self.assertEqual(len(labels), 1)
        self.assertTrue(8624790 in total_cluster_raw)
        self.assertTrue(3150921 in total_cluster_pf)
        self.assertTrue(472638150 in total_yield)

    def test_get_per_lane_sample_dist_plot(self):
        lane_plots = \
            get_per_lane_sample_dist_plot(
                sample_data=self.merged_sample_data)
        self.assertTrue(1 in lane_plots.keys())

    def test_get_project_summary_html_table(self):
        project_summary_data = \
            get_project_summary_html_table(
                sample_data=self.merged_sample_data)
        self.assertTrue(isinstance(project_summary_data, str))
        self.assertTrue('<td>IGFQ001220</td>' in project_summary_data)
        self.assertTrue('<td>95</td>' in project_summary_data)

    def test_get_per_lane_demult_table_data(self):
        table_data = \
            get_per_lane_demult_table_data(
                sample_data=self.merged_sample_data)
        self.assertTrue(1 in table_data.keys())
        #self.assertTrue('IGF118633' in table_data.get(1))
        #print(table_data.get(1).get('rows'))
        self.assertTrue('rows' in table_data.get(1))
        self.assertEqual(len(table_data.get(1).get('rows')), 95)

    def test_get_flowcell_project_summary_plot_for_db(self):
        project_summary_plot = \
            get_flowcell_project_summary_plot_for_db(
                summary_data=self.sum_df,
                sample_data=self.merged_sample_data)
        self.assertTrue('rows' in project_summary_plot.keys())
        self.assertEqual(len(project_summary_plot.get('rows')), 1)

    def test_get_undetermined_plot(self):
        undetermined_plots = \
            get_undetermined_plot(
                undetermined_data=self.undetermined_data)
        self.assertTrue(1 in undetermined_plots)
        self.assertTrue('labels' in undetermined_plots.get(1))
        self.assertTrue("GGCAAGTT+CCCCCCCC" in undetermined_plots.get(1).get('labels'))

    def test_get_undetermined_table(self):
        undetermined_tables = \
            get_undetermined_table(
                undetermined_data=self.undetermined_data)
        self.assertTrue(1 in undetermined_tables)
        self.assertTrue('rows' in undetermined_tables.get(1))
        self.assertEqual(len(undetermined_tables.get(1).get('rows')), 20)


if __name__ == '__main__':
  unittest.main()