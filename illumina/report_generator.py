import json
import os
import logging
import tempfile
from typing import Tuple, Optional
import pandas as pd
from datetime import datetime
from jinja2 import Environment, FileSystemLoader, select_autoescape
from iplotter import ChartJSPlotter
from iplotter import GCPlotter

def read_bcl2fastq_stats_data_from_pandas(data: dict) ->  Tuple[list, list, list]:
    '''
    A function for parsing Stats.json files from Illumina BCL2Fastq output

    :param data: A dictionary containing the following keys
        * ConversionResults
        * UnknownBarcodes
    :returns: Three lists
    '''
    try:
        row_l = list()
        row_s = list()
        unknown_df = list()
        for i in data.get('ConversionResults'):
            lane_number = i.get('LaneNumber')
            total_cluster_raw = i.get('TotalClustersRaw')
            total_cluster_pf = i.get('TotalClustersPF')
            total_yield = i.get('Yield')
            row_l.append({
                'Lane': lane_number,
                'Total_cluster_raw': total_cluster_raw,
                'Total_cluster_pf': total_cluster_pf,
                'Total_yield': total_yield})
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
        for unknown_entry in data.get('UnknownBarcodes'):
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
        raise ValueError(e)


def read_data_via_pandas(data_path: list) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    '''
    A function for reading list of Stats.json files from Illumina BCL2FASTQ output

    :param data_path: A list of Stats.json file paths
    :returns: Three Pandas DataFrames
    '''
    try:
        summary_records = pd.DataFrame()
        sample_records = pd.DataFrame()
        undetermined_records = pd.DataFrame()
        for f in data_path:
            with open(f, 'r') as fp:
                json_data = json.load(fp)
                row_l, row_s, unknown_df = \
                    read_bcl2fastq_stats_data_from_pandas(json_data)
                summary_records = \
                    pd.concat([summary_records, pd.DataFrame(row_l)], ignore_index=True)
                sample_records = \
                    pd.concat([sample_records, pd.DataFrame(row_s)], ignore_index=True)
                undetermined_records = \
                    pd.concat([undetermined_records, pd.DataFrame(unknown_df)], ignore_index=True)
        summary_records = \
            summary_records.groupby('Lane').agg(sum).\
            reset_index()
        return summary_records, sample_records, undetermined_records
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_stats_summary_table(
    sum_df: pd.DataFrame,
    lane_sample_df: pd.DataFrame) -> \
    pd.DataFrame:
    '''
    A function for calculating the summary table for de-multiplexing

    :param sum_df: A Pandas DataFrame containing summary data
    :param lane_sample_df: A Pandas DataFrame containing sample data
    :returns: A Pandas DataFrame containing the merged data table
    '''
    try:
        records = list()
        for (lane_id, sample_id), l_data in lane_sample_df.groupby(['Lane', 'Sample_ID']):
            sample_name = l_data['Sample_Name'].values[0]
            index_seq = l_data['Index_seq'].values[0]
            total_reads = l_data['Num_reads'].sum()
            perfect_barcodes = l_data['Perfect_barcode'].sum()
            total_yield = l_data['Yield'].sum()
            yield_q30 = l_data['Yield_q30'].sum()
            qual_score_sum = l_data['Qual_score_sum'].sum()
            total_pf_count = sum_df[sum_df['Lane'] == lane_id]['Total_cluster_pf'].values[0]
            records.append({
                'Lane': lane_id,
                'Sample_ID': sample_id,
                'Sample_Name': sample_name,
                'Barcode sequence': index_seq,
                'PF Clusters': '{:,}'.format(total_reads),
                '% of the lane': '{0:.2f}'.format(total_reads / total_pf_count * 100),
                '% Perfect barcode': '{0:.2f}'.format(perfect_barcodes / total_reads * 100),
                'Yield (Mbases)': '{:,.2f}'.format(total_yield / 10**6),
                '% >= Q30 bases': '{0:.2f}'.format(yield_q30 / total_yield * 100),
                'Mean Quality Score': '{0:.2f}'.format(qual_score_sum / total_yield)})
        records = pd.DataFrame(records)
        return records
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_samplesheet_records(samplesheets: list) -> pd.DataFrame:
    '''
    A function for parsing a list of samplesheet files

    :param samplesheets: A list of samplesheet files
    :returns: A Pandas DataFrame containing all the samplesheet data
    '''
    try:
        all_samplesheet_data = pd.DataFrame()
        for f in samplesheets:
            data_section = False
            samplesheet_data_list = list()
            with open(f, 'r') as fp:
                for i in fp:
                    if i.startswith('[Data]'):
                        data_section = True
                        continue
                    if data_section:
                        samplesheet_data_list.append(i.strip().split(','))
                samplesheet = \
                    pd.DataFrame(
                        samplesheet_data_list[1:],
                        columns=samplesheet_data_list[0])
                all_samplesheet_data = \
                    pd.concat(
                        [all_samplesheet_data, samplesheet],
                        ignore_index=True)
        return all_samplesheet_data
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_flowcell_summary_plots(
        summary_data: pd.DataFrame,
        div1_id: Optional[str] = 'chart1',
        div2_id: Optional[str] = 'chart2',
        width: Optional[int] = 400,
        height: Optional[int] = 400) -> \
        Tuple[str, str]:
    '''
    A function for plotting de-multiplexing summary

    :param summary_data: A Pandas DataFrame containing summary data
    :param div1_id: Div id for plot1, default 'chart1'
    :param div2_id: Div id for plot2, default 'chart2'
    :param width: Plot width, default 400
    :param height: Plot height, default 400
    :returns: Two strings containing the HTML plots
    '''
    try:
        summary_data.drop_duplicates(inplace=True)
        summary_data['Lane'] = summary_data['Lane'].astype(str)
        summary_data['Total_cluster_raw'] = summary_data['Total_cluster_raw'].astype(int)
        summary_data['Total_cluster_pf'] = summary_data['Total_cluster_pf'].astype(int)
        summary_data["Total_cluster_raw"] = summary_data["Total_cluster_raw"] /2
        ## generate flowcell summary plots
        data1 = {
            "labels": summary_data['Lane'].map(lambda x: 'Lane {0}'.format(x)).values.tolist(),
            "datasets": [{
                "label": "Total cluster raw",
                "data": summary_data["Total_cluster_raw"].astype(int).values.tolist(),
                "backgroundColor": 'rgba(255, 99, 132, 0.4)',
                "borderColor": 'rgba(255, 99, 132, 0.8)',
                "borderWidth": 1
            }, {
                "label": "Total cluster pf",
                "data": summary_data["Total_cluster_pf"].astype(int).values.tolist(),
                "backgroundColor": 'rgba(54, 162, 235, 0.4)',
                "borderColor": 'rgba(54, 162, 235, 0.8)',
                "borderWidth": 1
            }]
        }
        options1 = {
            "scales": {
                "y": {
                    "beginAtZero": "true"
                }
            },
            "responsive": "true",
            "title": {
                "display": "true",
                "text": 'Raw vs PF cluster counts per lane',
                "position": "bottom"
            },
            "plugins": {
                "legend": {
                    "position": 'top',
                }
            }
        }
        cj_plotter = ChartJSPlotter()
        summary_plot1 = \
            cj_plotter.render(
                data1,
                'bar',
                div_id=div1_id,
                options=options1,
                w=width,
                h=height)
        data2 = {
            "labels": summary_data['Lane'].map(lambda x: 'Lane {0}'.format(x)).values.tolist(),
            "datasets": [{
                "label": "Total yield",
                "data": summary_data["Total_yield"].astype(int).values.tolist(),
                "backgroundColor": 'rgba(75, 192, 192, 0.2)',
                "borderColor": 'rgba(75, 192, 192, 1)',
                "borderWidth": 1
            }]
        }
        options2 = {
            "scales": {
                "y": {
                    "beginAtZero": "true"
                }
            },
            "responsive": "true",
            "title": {
                "display": "true",
                "text": 'Total yield per lane',
                "position": "bottom"
            },
            "plugins": {
                "legend": {
                    "position": 'top',
                }
            }
        }
        summary_plot2 = \
            cj_plotter.render(
                data2,
                'bar',
                div_id=div2_id,
                options=options2,
                w=width,
                h=height)
        return summary_plot1, summary_plot2
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_per_lane_sample_dist_plot(
        sample_data: pd.DataFrame,
        bg_colors: list,
        border_colors: list,
        div_id_prefix: Optional[str] = 'chart_lane',
        plot_width: Optional[int] = 800,
        plot_height: Optional[int] = 400) -> \
        Tuple[dict, int]:
    '''
    A function for generating per lane sample distribution plot

    :param sample_data: A Pandas DataFrame containing the sample data
    :param bg_colors: A list of background colors
    :param border_colors: A list of border colors
    :param div_id_prefix: Div id prefix for the plots, default 'chart_lane'
    :param plot_width: Plot width, default 800
    :param plot_height: Plot height, default 400
    :returns: A dictionary containing the plot data and an integer for the recalculated plot height value
    '''
    try:
        # generate sample dis plots
        options = {
            "scales": {
                "y": {
                    "beginAtZero": "true"
                }
            },
            "responsive": "true",
            "plugins": {
                "legend": {
                    "position": 'top'}
            }}
        lane_plots = dict()
        max_samples = sample_data.groupby('Lane').agg(len)['Sample_ID'].max()
        if max_samples < 25:
            plot_height = 400
        elif max_samples >= 25 and max_samples < 50:
            plot_height = 800
        elif max_samples >= 50 and max_samples < 96:
            plot_height = 1000
        elif max_samples >= 96:
            plot_height = 1400
        for lane_id, l_data in sample_data.groupby('Lane'):
            lane_samples = l_data['Sample_ID'].values.tolist()
            datasets = list()
            counter = 0
            for project_id, p_data in l_data.groupby('Sample_Project'):
                pf_counts = \
                    p_data.\
                    set_index('Sample_ID')['PF Clusters'].\
                    reindex(lane_samples).\
                    fillna(0).\
                    values.tolist()
                datasets.append({
                    "label": project_id,
                    "data": pf_counts,
                    "backgroundColor": bg_colors[counter],
                    "borderColor": border_colors[counter],
                    "borderWidth": 1})
                counter += 1
            data = {
                "labels": lane_samples,
                "datasets": datasets}
            cj_plotter = ChartJSPlotter()
            plot_data = \
                cj_plotter.render(
                    data,
                    'horizontalBar',
                    div_id='{0}_{1}'.format(div_id_prefix, lane_id),
                    options=options,
                    w=plot_width,
                    h=plot_height)
            lane_plots.update({lane_id: plot_data})
        return lane_plots, plot_height
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_project_summary_html_table(sample_data: pd.DataFrame) -> str:
    try:
        project_summary_data = \
            sample_data.\
                groupby(['Lane', 'Sample_Project']).\
                agg(len)['Sample_ID'].reset_index().\
                to_html(
                    index=False,
                    classes=['table', 'table-sm', 'table-hover'],
                    justify='center',
                    border=0)
        return project_summary_data
    except Exception as e:
        raise ValueError(e)


def get_per_lane_demult_table_data(sample_data: pd.DataFrame) -> dict:
    try:
        data = list()
        data.append([
            'Lane',
            'Sample_ID',
            'Sample_Name',
            'Sample_Project',
            'Barcode sequence',
            'PF Clusters',
            '% of the lane',
            '% Perfect barcode',
            'Yield (Mbases)',
            '% >= Q30 bases',
            'Mean Quality Score'])
        data.extend(sample_data[[
            'Lane',
            'Sample_ID',
            'Sample_Name',
            'Sample_Project',
            'Barcode sequence',
            'PF Clusters',
            '% of the lane',
            '% Perfect barcode',
            'Yield (Mbases)',
            '% >= Q30 bases',
            'Mean Quality Score']].\
            values.tolist())
        table_data = dict()
        for lane_id, l_data in pd.DataFrame(data[1:], columns=data[0]).groupby('Lane'):
            l_data = l_data.fillna('')
            l_json_data = [l_data.columns.tolist()]
            l_json_data.extend(l_data.values.tolist())
            table_data.update({lane_id: json.dumps(l_json_data)})
        return table_data
    except Exception as e:
        raise ValueError(e)


def get_flowcell_project_summary_plot(
        summary_data: pd.DataFrame,
        sample_data: pd.DataFrame,
        div_id: Optional[str] = 'project_summary_plot',
        width: Optional[int] = 900,
        height: Optional[int] = 700) -> str:
    try:
        lane_project_pf_counts = \
            sample_data.\
                groupby(['Lane', 'Sample_Project']).\
                agg(sum).\
                reset_index()[['Lane', 'Sample_Project', 'PF Clusters']]
        summary_data['Lane'] = summary_data['Lane'].astype(int)
        lane_project_pf_counts['Lane'] = \
            lane_project_pf_counts['Lane'].astype(int)
        m_df = \
            lane_project_pf_counts.\
                set_index('Lane').\
                join(summary_data.set_index('Lane'), how='left')
        ld = list()
        for lane_id, l_data in m_df.groupby('Lane'):
            lane_pf_count = 0 
            row_data = {'Lane': 'Lane {0}'.format(lane_id)}
            for project_id, p_data in l_data.groupby('Sample_Project'):
                lane_pf_count += p_data['PF Clusters'].values[0]
                row_data.update({project_id: p_data['PF Clusters'].sum()})
                row_data.update({'undetermined': l_data['Total_cluster_pf'].values[0] - lane_pf_count})
            ld.append(row_data)
        data5 = list()
        data5.append(pd.DataFrame(ld).columns.tolist())
        data5.extend(pd.DataFrame(ld).values.tolist())
        options5 = {
            "title": 'Project summary plot',
            "width": width,
            "height": height,
            "chartArea": {"left": 150, "width": "60%"},
            "legend": {"position": 'right', "maxLines": 20, "fontSize": 5},
            "dataOpacity": 0.5,
            "colors": [
                '#0173B2', '#DE8F05', '#029E73', '#D55E00', '#CC78BC', '#CA9161',
                '#FBAFE4', '#949494', '#ECE133', '#56B4E9', '#0173B2'],
            "bar": {"groupWidth": '70%'},
            "isStacked": "percent"}
        gcplotter = GCPlotter()
        project_summary_plot = \
            gcplotter.render(
                data5,
                chart_type="BarChart",
                div_id=div_id,
                chart_package='corechart',
                options=options5)
        return project_summary_plot
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_undetermined_plot(
        undetermined_data: pd.DataFrame,
        bg_colors: list,
        border_colors: list,
        div_id_prefix: Optional[str] = 'ud_chart_lane',
        barcode_limit: Optional[int] = 20,
        width: Optional[int] = 600,
        height: Optional[int] = 200) -> dict:
    try:
        # generating undetermined plots
        options6 = {
            "scales": {
                "y": {
                    "beginAtZero": "true"
                }
            },
            "responsive": "true",
            "plugins": {
                "legend": {
                    "position": 'top'}
            }
        }
        udf = pd.DataFrame()
        counter = 0
        undetermined_plots = dict()
        for lane_id, l_data in undetermined_data.groupby('Lane'):
            sorted_df = l_data.sort_values('Reads', ascending=False).head(barcode_limit).copy()
            if len(udf.index) > 0:
                udf = pd.concat([udf, sorted_df], ignore_index=True)
            else:
                udf = sorted_df.copy()
            barcode_labels = sorted_df['Barcode'].values.tolist()
            barcode_count = sorted_df['Reads'].values.tolist()
            data6 = {
                "labels": barcode_labels,
                "datasets": [{
                    "label": "Lane {0}".format(lane_id),
                    "data": barcode_count,
                    "backgroundColor": bg_colors[counter],
                    "borderColor": border_colors[counter],
                    "borderWidth": 1}]
            }
            counter += 1
            cj_plotter = ChartJSPlotter()
            plot = \
                cj_plotter.render(
                    data6,
                    'horizontalBar',
                    div_id='{0}_{1}'.format(div_id_prefix, lane_id),
                    options=options6,
                    w=width,
                    h=height)
            undetermined_plots.update({lane_id: plot})
        return undetermined_plots
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def get_undetermined_table(undetermined_data: pd.DataFrame) -> dict:
    try:
        udf = pd.DataFrame()
        undetermined_data['Reads'] = \
            undetermined_data['Reads'].astype(int)
        for lane_id, l_data in undetermined_data.groupby('Lane'):
            sorted_df = l_data.sort_values('Reads', ascending=False).head(20).copy()
            if len(udf.index) > 0:
                udf = pd.concat([udf, sorted_df], ignore_index=True)
            else:
                udf = sorted_df.copy()
        udf['Barcode_I2_RC'] = \
            udf['Barcode'].\
            map(lambda x: \
                '+'.join([
                    x.split('+')[0], x.split('+')[1].translate(str.maketrans('ATGC','TACG'))[::-1]]) \
                        if '+' in x else '')
        udf['Barcode_I1_RC'] = \
            udf['Barcode'].\
                map(lambda x: \
                    '+'.join([
                        x.split('+')[0].translate(str.maketrans('ATGC','TACG'))[::-1], x.split('+')[1]]) \
                            if '+' in x else x.translate(str.maketrans('ATGC','TACG'))[::-1])
        u_data = list()
        u_data.append([
            'Lane',
            'Barcode',
            'Barcode_I1_RC',
            'Barcode_I2_RC',
            'Reads'])
        u_data.extend(udf[[
            'Lane',
            'Barcode',
            'Barcode_I1_RC',
            'Barcode_I2_RC',
            'Reads']].\
            values.tolist())
        undetermined_tables = dict()
        for lane_id, l_data in pd.DataFrame(u_data[1:], columns=u_data[0]).groupby('Lane'):
            l_data = l_data.fillna('')
            l_json_data = [l_data.columns.tolist()]
            l_json_data.extend(l_data.values.tolist())
            undetermined_tables.update({lane_id: json.dumps(l_json_data)})
        return undetermined_tables
    except Exception as e:
        logging.error(e)
        raise ValueError(e)


def create_demux_html_report(
        project_summary_plot: str,
        project_summary_data: str,
        summary_plot1: str,
        summary_plot2: str,
        plot_height: int,
        lane_plots: list,
        sample_tables: dict,
        undetermined_plots: list,
        undetermined_tables: dict,
        template_path: str,
        output_file: str,
        seqrun_id: str) -> None:
    try:
        template_env = \
            Environment(
                loader=FileSystemLoader(searchpath=os.path.dirname(template_path)),
                autoescape=select_autoescape(['html', 'xml']))
        html = \
            template_env.\
            get_template(os.path.basename(template_path))
        html.\
            stream(dict(
                DATE_STAMP=datetime.strftime(datetime.now(),'%Y %h %d - %H:%M'),
                SEQRUN_ID=seqrun_id,
                PROJECT_SUMMARY_PLOT=project_summary_plot,
                PROJECT_SUMMARY=project_summary_data,
                SUMMARY_PLOT1=summary_plot1,
                SUMMARY_PLOT2=summary_plot2,
                PLOT_HEIGHT=plot_height,
                LANE_PLOTS=lane_plots,
                SAMPLE_TABLES=sample_tables,
                UNDETERMINED_PLOTS=undetermined_plots,
                UNDETERMINED_TABLES=undetermined_tables)).\
            dump(output_file)
    except Exception as e:
        logging.error(e)
        raise ValueError(e)

def combine_data_and_create_report(
    sum_df: pd.DataFrame,
    lane_sample_df: pd.DataFrame,
    undetermined_data: pd.DataFrame,
    samplesheets: list,
    seqrun_id: str,
    template_path: str,
    output_file: str) -> None:
    try:
        bg_colors = [
            'rgba(255, 99, 132, 0.4)',
            'rgba(54, 162, 235, 0.4)',
            'rgba(255, 206, 86, 0.4)',
            'rgba(75, 192, 192, 0.4)',
            'rgba(153, 102, 255, 0.4)',
            'rgba(255, 159, 64, 0.4)',
            'rgba(255, 159, 10, 0.4)',
            'rgba(255, 159, 192, 0.4)']
        border_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        all_samplesheet_data = \
            get_samplesheet_records(samplesheets)
        sample_data = \
            get_stats_summary_table(
                sum_df,
                lane_sample_df)
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
        summary_plot1, summary_plot2 = \
            get_flowcell_summary_plots(
                summary_data=sum_df)
        lane_plots, plot_height = \
            get_per_lane_sample_dist_plot(
                sample_data=merged_sample_data,
                bg_colors=bg_colors,
                border_colors=border_colors)
        project_summary_data = \
            get_project_summary_html_table(
                sample_data=merged_sample_data)
        table_data = \
            get_per_lane_demult_table_data(
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


def prepare_report_using_pandas(
        data_path: str,
        samplesheets: str,
        seqrun_id: str,
        template_path: str,
        output_file: str) -> None:
    try:
        if not os.path.exists(data_path) or \
           not os.path.exists(samplesheets):
           raise IOError('Stats.josn or SampleSheet.csv files not found')
        with open(data_path, 'r') as fp:
            stats_jsons = [f.strip() for f in fp]
        with open(samplesheets, 'r') as fp:
            samplesheet_csvs = [f.strip() for f in fp]
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_pandas(
                data_path=stats_jsons)
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