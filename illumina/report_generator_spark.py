import logging, os
import pandas as pd
from typing import Tuple, Optional
from pyspark.sql import SparkSession
from pyspark.sql.types import StructField, StructType, StringType, IntegerType
from pyspark.sql.types import ArrayType, LongType, BooleanType, MapType
from pyspark.sql.functions import explode, col, split, regexp_replace, sum
from .report_generator import combine_data_and_create_report

def get_local_spark_session(
        threads: Optional[int]=1,
        executor_memory_in_gb: Optional[int]=4) -> SparkSession:
    try:
        spark = \
            SparkSession.\
            builder.\
            config("spark.executor.cores", threads).\
            config("spark.sql.execution.arrow.pyspark.enabled", "true").\
            config("spark.executor.memory", "{0}G".format(executor_memory_in_gb)).\
            getOrCreate()
        return spark
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_json_schema() -> StructType():
    try:
        schema = \
            StructType([
                StructField("Flowcell", StringType(), True),
                StructField("RunNumber", IntegerType(), True),
                StructField("RunId", StringType(), True),
                StructField("ReadInfosForLanes",
                    ArrayType(
                        StructType([
                            StructField("LaneNumber", IntegerType(), True),
                            StructField("ReadInfos",
                                ArrayType(
                                    StructType([
                                        StructField("Number", IntegerType(), True),
                                        StructField("NumCycles", IntegerType(), True),
                                        StructField("IsIndexedRead", BooleanType(), True)])),
                                True)])),
                    True),
                StructField("ConversionResults",
                    ArrayType(
                        StructType([
                            StructField("LaneNumber", IntegerType(), True),
                            StructField("TotalClustersRaw", LongType(), True),
                            StructField("TotalClustersPF", LongType(), True),
                            StructField("Yield", LongType(), True),
                            StructField("DemuxResults",
                                ArrayType(
                                    StructType([
                                        StructField("SampleId", StringType(), True),
                                        StructField("SampleName", StringType(), True),
                                        StructField("IndexMetrics",
                                            ArrayType(
                                                StructType([
                                                    StructField("IndexSequence", StringType(), True),
                                                    StructField("MismatchCounts",
                                                        MapType(StringType(), IntegerType(), True),
                                                    True)
                                            ])),
                                        True),
                                        StructField("NumberReads", IntegerType(), True),
                                        StructField("Yield", LongType(), True),
                                        StructField("ReadMetrics",
                                            ArrayType(
                                                StructType([
                                                    StructField("ReadNumber", IntegerType(), True),
                                                    StructField("Yield", LongType(), True),
                                                    StructField("YieldQ30", LongType(), True),
                                                    StructField("QualityScoreSum", LongType(), True),
                                                    StructField("TrimmedBases", IntegerType(), True)
                                            ])),
                                            True)
                                ])),
                                True),
                            StructField("Undetermined",
                                StructType([
                                    StructField("NumberReads", IntegerType(), True),
                                    StructField("Yield", LongType(), True),
                                    StructField("ReadMetrics",
                                        ArrayType(
                                            StructType([
                                                StructField("ReadNumber", IntegerType(), True),
                                                StructField("Yield", LongType(), True),
                                                StructField("YieldQ30", LongType(), True),
                                                StructField("QualityScoreSum", LongType(), True),
                                                StructField("TrimmedBases", IntegerType(), True)
                                        ])),
                                        True)
                                ]),
                                True),
                    ])),
                    True),
                StructField("UnknownBarcodes",
                    ArrayType(
                        MapType(StringType(), StringType(), True)),
                    True)
            ])
        return schema
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def read_data_via_pyspark(
        spark: SparkSession,
        data_files: list) -> \
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    try:
        schema = get_json_schema()
        df = \
            spark.\
            read.\
            format("json").\
            option("mode", "failfast").\
            option('inferSchema', 'false').\
            option('multiLine', 'true').\
            schema(schema).\
            load(data_files)
        sum_df = \
            df.\
            withColumn('ConversionResultsExploded', explode('ConversionResults')).\
            select(
                'RunId',
                'ConversionResultsExploded.LaneNumber',
                'ConversionResultsExploded.TotalClustersRaw',
                'ConversionResultsExploded.TotalClustersPF',
                'ConversionResultsExploded.Yield').\
            groupby(['RunId', 'LaneNumber']).\
            agg(
                sum('TotalClustersRaw').alias('TotalClustersRaw'),
                sum('TotalClustersPF').alias('TotalClustersPF'),
                sum('Yield').alias('Yield')).\
            select(
                col('LaneNumber').alias('Lane'),
                col('TotalClustersRaw').cast('long').alias('Total_cluster_raw'),
                col('TotalClustersPF').cast('long').alias('Total_cluster_pf'),
                col('Yield').cast('long').alias('Total_yield')).\
            toPandas()
        lane_sample_df = \
            df.\
            withColumn('ConversionResultsExploded', explode('ConversionResults')).\
            withColumn('DemuxResultsExploded', explode('ConversionResults.DemuxResults')).\
            withColumn('DemuxResultsExplodedRe', explode('DemuxResultsExploded')).\
            withColumn('IndexMetricsExploded', explode('DemuxResultsExplodedRe.IndexMetrics')).\
            withColumn('ReadMetricsExploded', explode('DemuxResultsExplodedRe.ReadMetrics')).\
            selectExpr(
                'RunId',
                'ConversionResultsExploded.LaneNumber as Lane',
                'IndexMetricsExploded.IndexSequence as Index_seq',
                'IndexMetricsExploded.MismatchCounts[0] as Perfect_barcode',
                'IndexMetricsExploded.MismatchCounts[1] as OneMismatchBarcode',
                'DemuxResultsExplodedRe.SampleId as Sample_ID',
                'DemuxResultsExplodedRe.SampleName as Sample_Name',
                'DemuxResultsExplodedRe.NumberReads as PFClusters',
                'ReadMetricsExploded.ReadNumber as Num_reads',
                'ReadMetricsExploded.Yield as Yield',
                'ReadMetricsExploded.YieldQ30 as YieldQ30',
                'ReadMetricsExploded.QualityScoreSum as QualityScoreSum',
                'ReadMetricsExploded.TrimmedBases as TrimmedBases').\
            groupby([
                'RunId',
                'Lane',
                'Index_seq',
                'Sample_ID',
                'Sample_Name',
                'PFClusters',
                'Perfect_barcode']).\
            agg(
                sum('Yield').alias('Yield'),
                sum('YieldQ30').alias('Yield_q30'),
                sum('QualityScoreSum').alias('Qual_score_sum')).\
            select(
                'Lane',
                'Sample_ID',
                'Sample_Name',
                col('PFClusters').alias('Num_reads'),
                'Index_seq',
                'Perfect_barcode',
                'Yield',
                'Yield_q30',
                'Qual_score_sum').\
            toPandas()
        undetermined_data = \
            df.\
            withColumn('UnknownBarcodesExploded', explode('UnknownBarcodes')).\
            select(
                'UnknownBarcodesExploded.Lane',
                'UnknownBarcodesExploded.Barcodes').\
            select(
                'Lane',
                regexp_replace(col('Barcodes'), "\{|\}|\"", "").alias('Barcodes')).\
            withColumn('sb', split('Barcodes', ',')).\
            withColumn('sb_1', explode('sb')).\
            select(
                'Lane',
                split('sb_1', ':').getItem(0).alias('Barcode'),
                split('sb_1', ':').getItem(1).alias('Reads')).\
            toPandas()
        return sum_df, lane_sample_df, undetermined_data
    except Exception as e:
        raise ValueError(e)

def prepare_report_using_spark(
        data_path: str,
        samplesheets: str,
        seqrun_id: str,
        template_path: str,
        output_file: str,
        threads: Optional[int] = 4,
        executor_memory_in_gb: Optional[int] = 4) -> None:
    try:
        if not os.path.exists(data_path) or \
           not os.path.exists(samplesheets):
           raise IOError('Stats.josn or SampleSheet.csv files not found')
        spark = \
            get_local_spark_session(
                threads=threads,
                executor_memory_in_gb=executor_memory_in_gb)
        with open(data_path, 'r') as fp:
            stats_jsons = [f.strip() for f in fp]
        with open(samplesheets, 'r') as fp:
            samplesheet_csvs = [f.strip() for f in fp]
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_pyspark(
                spark=spark,
                data_files=stats_jsons)
        combine_data_and_create_report(
            sum_df=sum_df,
            lane_sample_df=lane_sample_df,
            undetermined_data=undetermined_data,
            samplesheets=samplesheet_csvs,
            seqrun_id=seqrun_id,
            template_path=template_path,
            output_file=output_file)
    except Exception as e:
        logging.error(e)
        raise ValueError(e)