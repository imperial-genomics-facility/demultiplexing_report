import logging
from typing import Tuple, Optional
import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client
from dask.distributed import LocalCluster
from .report_generator import combine_data_and_create_report


def read_bcl2fastq_stats_data_using_dask(x: dict) -> Tuple[list, list, list]:
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


def get_local_dask_cluster(
        n_workers: Optional[int] = 4,
        threads_per_worker: Optional[int] = 1,
        memory_limit: Optional[str] = '1GB') -> LocalCluster:
    try:
        cluster = \
            LocalCluster(
                n_workers=n_workers,
                threads_per_worker=threads_per_worker,
                memory_limit=memory_limit)
        return cluster
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def read_data_via_dask_cluster(
        data_path: list,
        dask_cluster: LocalCluster) -> \
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    try:
        client = \
            Client(dask_cluster)
        df_all = dd.read_json(data_path, orient='columns')
        df_all = df_all.compute()
        all_pdf = \
            df_all.apply(
                read_bcl2fastq_stats_data_using_dask,
                axis=1,
                result_type='expand')
        lane_df = \
            pd.DataFrame([v[0] for v in all_pdf[0].values])
        sum_df = \
            lane_df.groupby('Lane').agg(sum)
        sum_df = \
            sum_df.reset_index()
        lane_sample_df = \
            pd.DataFrame([w for v in all_pdf[1].values for w in v])
        lane_unknown_df = \
            pd.DataFrame([w for v in all_pdf[2].values for w in v])
        undetermined_data = \
            lane_unknown_df.groupby(['Lane','Barcode']).agg(sum)
        undetermined_data = \
            undetermined_data.reset_index()
        client.close()
        return sum_df, lane_sample_df, undetermined_data
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def prepare_report_using_dask(
        data_path: list,
        samplesheets: list,
        seqrun_id: str,
        template_path: str,
        output_file: str,
        n_workers: Optional[int] = 4,
        threads_per_worker: Optional[int] = 1,
        memory_limit: Optional[str] = '1GB') -> None:
    try:
        cluster = \
            get_local_dask_cluster(
                n_workers=n_workers,
                threads_per_worker=threads_per_worker,
                memory_limit=memory_limit)                                      # get dask cluster
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_dask_cluster(
                data_path=data_path)                                            # read files using dask cluster
        cluster.close()                                                         # close dask cluster
        combine_data_and_create_report(
            sum_df=sum_df,
            lane_sample_df=lane_sample_df,
            undetermined_data=undetermined_data,
            samplesheets=samplesheets,
            seqrun_id=seqrun_id,
            template_path=template_path,
            output_file=output_file)
    except Exception as e:
        logging.error(e)
        raise ValueError(e)