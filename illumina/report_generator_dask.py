import logging
import typing
from typing import List
import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client
from report_generator import get_samplesheet_records
from report_generator import get_stats_summary_table
from report_generator import get_flowcell_summary_plots
from report_generator import get_per_lane_sample_dist_plot
from report_generator import get_project_summary_html_table
from report_generator import get_demult_per_lane_demult_table_data
from report_generator import get_flowcell_project_summary_plot
from report_generator import get_undetermined_plot
from report_generator import get_undetermined_table
from report_generator import create_demux_html_report

def read_bcl2fastq_stats_data_using_dask(x: dict) -> List[list, list, list]:
    try:
        i = x.get('ConversionResults')
        lane_number = i.get('LaneNumber')
        total_cluster_raw = i.get('TotalClustersRaw')
        total_cluster_pf = i.get('TotalClustersPF')
        total_yield = i.get('Yield')
        row_l = [{
            'Lane': lane_number,
            'Total_cluster_raw': total_cluster_raw,
            'Total_cluster_pf': total_cluster_pf,
            'Total_yield': total_yield}]
        row_s = list()
        demux_results = i.get('DemuxResults')
        for j in demux_results:
            sample_id = j.get('SampleId')
            sample_name = j.get('SampleName')
            index = j.get('IndexMetrics')[0].get('IndexSequence')
            num_reads = j.get('NumberReads')
            yield_val = j.get('Yield')
            perfect_barcodes = j['IndexMetrics'][0]['MismatchCounts']['0']
            yield_q30 = 0
            qual_score_sum = 0
            read_metrics = j.get('ReadMetrics')
            for read in read_metrics:
                q30_bases = int(read.get('YieldQ30'))
                yield_q30 += q30_bases
                qual_score = int(read.get('QualityScoreSum'))
                qual_score_sum += qual_score
            row_s.append({
                'Lane': lane_number,
                'Sample_ID': sample_id,
                'Sample_Name': sample_name,
                'Index_seq': index,
                'Num_reads': num_reads,
                'Perfect_barcode': perfect_barcodes,
                'Yield_q30': yield_q30,
                'Yield': int(yield_val),
                'Qual_score_sum': qual_score_sum})
        unknown_df = list()
        unknown_entry = x['UnknownBarcodes']
        lane_id = unknown_entry.get('Lane')
        barcodes = unknown_entry.get('Barcodes')
        for barcode, read in barcodes.items():
            unknown_df.\
            append({
                'Lane': lane_id,
                'Barcode': barcode,
                'Reads': read })
        return row_l, row_s, unknown_df
    except Exception as e:
        logging.errror(e)
        raise ValueError(e)


def read_data_via_dask_cluster(
        data_path:list, n_workers:int=4, threads_per_worker:int=1, memory_limit:str='1GB') -> \
        List[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    try:
        client = \
            Client(
                n_workers=n_workers,
                threads_per_worker=threads_per_worker,
                memory_limit=memory_limit)
        df_all = dd.read_json(data_path, orient='columns')
        df_all = df_all.compute()
        all_pdf = df_all.apply(read_bcl2fastq_stats_data_using_dask, axis=1, result_type='expand')
        client.close()
        lane_df = pd.DataFrame([v[0] for v in all_pdf[0].values])
        sum_df = lane_df.groupby('Lane').agg(sum)
        sum_df = sum_df.reset_index()
        lane_sample_df = pd.DataFrame([w for v in all_pdf[1].values for w in v])
        lane_unknown_df = pd.DataFrame([w for v in all_pdf[2].values for w in v])
        undetermined_data = lane_unknown_df.groupby(['Lane','Barcode']).agg(sum)
        undetermined_data = undetermined_data.reset_index()
        return sum_df, lane_sample_df, undetermined_data
    except Exception as e:
        logging.error(e)
        raise ValueError(e)

def prepare_report_using_dask(
        data_path:list, samplesheets:list, seqrun_id:str, template_path:str, output_file:str) -> None:
    try:
        bg_colors = [
            'rgba(255, 99, 132, 0.4)',
            'rgba(54, 162, 235, 0.4)',
            'rgba(255, 206, 86, 0.4)',
            'rgba(75, 192, 192, 0.4)',
            'rgba(153, 102, 255, 0.4)',
            'rgba(255, 159, 64, 0.4)',
            'rgba(255, 159, 10, 0.4)',
            'rgba(255, 159, 192, 0.4)'
        ]
        border_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)'
        ]
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_dask_cluster(data_path=data_path)
        all_samplesheet_data = get_samplesheet_records(samplesheets)
        sample_data = get_stats_summary_table(sum_df, lane_sample_df)
        merged_sample_data = \
            sample_data.set_index('Sample_ID').\
            join(all_samplesheet_data.set_index('Sample_ID')['Sample_Project'], how='left').\
            reset_index().\
            drop_duplicates()
        merged_sample_data['PF Clusters'] = \
            merged_sample_data['PF Clusters'].map(lambda x: x.replace(',',''))
        merged_sample_data['PF Clusters'] = merged_sample_data['PF Clusters'].astype(int)
        summary_plot1, summary_plot2 = \
            get_flowcell_summary_plots(
                summary_data=sum_df,
                bg_colors=bg_colors,
                border_colors=border_colors)
        lane_plots, plot_height = \
            get_per_lane_sample_dist_plot(
                sample_data=merged_sample_data,
                bg_colors=bg_colors,
                border_colors=border_colors)
        project_summary_data = \
            get_project_summary_html_table(
                sample_data=merged_sample_data)
        table_data = \
            get_demult_per_lane_demult_table_data(
                sample_data=merged_sample_data)
        project_summary_plot = \
            get_flowcell_project_summary_plot(
                summary_data=sum_df,
                sample_data=merged_sample_data)
        undetermined_plots = \
            get_undetermined_plot(
                undetermined_data=undetermined_data,
                bg_colors=bg_colors,
                border_colors=border_colors)
        undetermined_tables = \
            get_undetermined_table(
                undetermined_data=undetermined_data)
        create_demux_html_report(
            project_summary_plot=project_summary_plot,
            project_summary_data=project_summary_data,
            summary_plot1=summary_plot1,
            summary_plot2=summary_plot2,
            plot_height=plot_height,
            lane_plots=lane_plots,
            sample_tables=table_data,
            undetermined_plots=undetermined_plots,
            undetermined_tables=undetermined_tables,
            template_path=template_path,
            output_file=output_file,
            seqrun_id=seqrun_id)
    except Exception as e:
        logging.error(e)
        raise ValueError(e)