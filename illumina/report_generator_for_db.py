import typing
import pandas as pd
from typing import Tuple
import os, json, gviz_api, tempfile
from shutil import copy2
from illumina.report_generator import (
    read_data_via_pandas,
    get_stats_summary_table,
    get_samplesheet_records)

def get_flowcell_summary_data(summary_data):
    """
    A function for flowcell summary table

    :param summary_data: A Pandas DataFrame containing the summary data
    :returns: Four lists containing plot label, raw cluster datata, PF cluster data and yield data
    """
    try:
        summary_data.drop_duplicates(inplace=True)
        summary_data['Lane'] = summary_data['Lane'].astype(str)
        summary_data['Total_cluster_raw'] = summary_data['Total_cluster_raw'].astype(int)
        summary_data['Total_cluster_pf'] = summary_data['Total_cluster_pf'].astype(int)
        summary_data["Total_cluster_raw"] = summary_data["Total_cluster_raw"] /2
        labels = \
            summary_data['Lane'].map(lambda x: 'Lane {0}'.format(x)).values.tolist()
        total_cluster_raw = \
            summary_data["Total_cluster_raw"].astype(int).values.tolist()
        total_cluster_pf = \
            summary_data["Total_cluster_pf"].astype(int).values.tolist()
        total_yield = \
            summary_data["Total_yield"].astype(int).values.tolist()
        return labels, total_cluster_raw, total_cluster_pf, total_yield
    except Exception as e:
        raise ValueError("Failed to get flowcell summary data, error: {0}".format(e))


def get_per_lane_sample_dist_plot(sample_data, bg_colors, border_colors):
    """
    A function for returing sample distribution plots

    :params sample_data: A Pandas DataFrame containing sample data
    :returns: A dictionary containing sample distribution plots
    """
    try:
        lane_plots = dict()
        for lane_id, l_data in sample_data.groupby('Lane'):
            lane_samples = l_data['Sample_ID'].values.tolist()
            datasets = list()
            counter = 0
            project_count = \
                len(l_data.groupby('Sample_Project').groups.keys())
            if project_count > len(bg_colors) or \
               project_count > len(border_colors):
                raise ValueError("Not enough color codes for project")
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
                    "borderColor": border_colors[counter]})
                counter += 1
            data = {
                "labels": lane_samples,
                "datasets": datasets}
            lane_plots.update({lane_id: data})
        return lane_plots
    except Exception as e:
        raise ValueError("Failed to get sample distribution data, error: {0}".format(e))


def get_project_summary_html_table(sample_data):
    """
    A function for generating project summary html table

    :param sample_data: A Pandas DataFrame containing sample data
    :returns: A string containing the html table data
    """
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
        raise ValueError("Failed to get flowcell summary table, error: {0}".format(e))


def get_per_lane_demult_table_data(sample_data):
    """
    A function for generating de-multiplexing table data

    :param sample_data: A Pandas DataFrame containing sample data
    :returns: A dictionary containing the gviz table data
    """
    try:
        data = list()
        columns = [
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
            'Mean Quality Score']
        data.append(columns)
        data.extend(sample_data[columns].\
            values.tolist())
        table_data = dict()
        description = dict()
        for entry in columns:
            description.update({entry: ("string", entry)})
        for lane_id, l_data in pd.DataFrame(data[1:], columns=data[0]).groupby('Lane'):
            l_data = l_data.fillna('')
            gviz_data = gviz_api.DataTable(description)
            gviz_data.LoadData(l_data.to_dict(orient="records"))
            table_data.update({lane_id: json.loads(gviz_data.ToJSon(columns_order=columns))})
        return table_data
    except Exception as e:
        raise ValueError("Failed to get de-multiplexing table data, error: {0}".format(e))


def get_flowcell_project_summary_plot_for_db(summary_data, sample_data):
    """
    A function for getting flowcell summary plots

    :param summary_data: A Pandas DataFrame containing flowcell summary data
    :param sample_data: A Pandas DataFrame containing samples data
    :returns: A dictionary containing the plot data
    """
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
        m_df.fillna(0, inplace=True)
        ld = list()
        for lane_id, l_data in m_df.groupby('Lane'):
            lane_pf_count = 0
            row_data = {'Lane': 'Lane {0}'.format(lane_id)}
            for project_id, p_data in l_data.groupby('Sample_Project'):
                lane_pf_count += p_data['PF Clusters'].values[0]
                row_data.update({project_id: p_data['PF Clusters'].sum()})
                row_data.update({'undetermined': l_data['Total_cluster_pf'].values[0] - lane_pf_count})
            ld.append(row_data)
        project_summary_plot = list()
        project_summary_plot.append(pd.DataFrame(ld).columns.tolist())
        project_summary_plot.extend(pd.DataFrame(ld).fillna(0).values.tolist())
        description = {
            c:("string",c) if c=='Lane' else ("number",c)
                for c in pd.DataFrame(ld).columns.tolist()}
        data_table = gviz_api.DataTable(description)
        data_table.\
            LoadData(pd.DataFrame(ld).fillna(0).to_dict(orient="records"))
        plot_data = \
            json.loads(data_table.ToJSon(columns_order=pd.DataFrame(ld).columns.tolist()))
        return plot_data
    except Exception as e:
        raise ValueError("Failed to get flowcell summary plot, error: {0}".format(e))


def get_undetermined_table(undetermined_data):
    """
    A function for creating undetermined data table

    :param undetermined_data: A Pandas DataFrame containing the undetermined data
    :returns: A dictionary containing the html table data for individual lanes
    """
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
        columns = [
            'Lane',
            'Barcode',
            'Barcode_I1_RC',
            'Barcode_I2_RC',
            'Reads']
        u_data.append(columns)
        u_data.extend(udf[columns].\
            values.tolist())
        undetermined_tables = dict()
        description = {c: ("string", c) for c in columns}
        for lane_id, l_data in pd.DataFrame(u_data[1:], columns=u_data[0]).groupby('Lane'):
            l_data = l_data.fillna('')
            gviz_data = gviz_api.DataTable(description)
            gviz_data.LoadData(l_data.to_dict(orient="records"))
            undetermined_tables.update({lane_id: json.loads(gviz_data.ToJSon(columns_order=columns))})
        return undetermined_tables
    except Exception as e:
        raise ValueError("Failed to get undetermined table data, error: {0}".format(e))


def get_undetermined_plot(undetermined_data, barcode_limit=20):
    """
    A function got getting undteremined plot data

    :param undetermined_data: A Pandas DataFrame containing undetermined data
    :param barcode_limit: Barcode limit, default 20
    :returns: A dictionary containing the plot data for individual lanes
    """
    try:
        udf = pd.DataFrame()
        undetermined_plots = dict()
        for lane_id, l_data in undetermined_data.groupby('Lane'):
            sorted_df = l_data.sort_values('Reads', ascending=False).head(barcode_limit).copy()
            if len(udf.index) > 0:
                udf = pd.concat([udf, sorted_df], ignore_index=True)
            else:
                udf = sorted_df.copy()
            barcode_labels = sorted_df['Barcode'].values.tolist()
            barcode_count = sorted_df['Reads'].values.tolist()
            data = {
                "labels": barcode_labels,
                "data": barcode_count}
            undetermined_plots.update({lane_id: data})
        return undetermined_plots
    except Exception as e:
        raise ValueError("Failed to get undetermined plot data, error: {0}".format(e))


def create_plot_json_for_database(
        run_name, samplesheet_tag, stat_files,
        samplesheet_files, output_dir):
    """
    A function for creating json de-multiplexing report data

    :param run_name: A string containing the run_name
    :param samplesheet_tag: A string containing the samplesheet tag
    :param stat_files: A file containing the de-multiplexing Stats.json files
    :param samplesheet_files: A list of samplesheet files that used for de-multiplexing
    :param output_dir: A string containing the output dir
    :returns: None
    """
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
        final_json_file = \
            os.path.join(
                output_dir,
                '{0}_{1}.json'.format(run_name, samplesheet_tag))
        os.makedirs(output_dir, exist_ok=True)
        if os.path.exists(final_json_file):
            raise IOError(
                    'Output file {0} already present. Remove it before re-run.'.\
                        format(final_json_file))
        with open(stat_files, 'r') as fp:
            stats_jsons = [f.strip() for f in fp]
        with open(samplesheet_files, 'r') as fp:
            samplesheet_list = [f.strip() for f in fp]
        sum_df, lane_sample_df, undetermined_data = \
            read_data_via_pandas(
                data_path=stats_jsons)
        sample_data = \
            get_stats_summary_table(
                sum_df,
                lane_sample_df)
        all_samplesheet_data = \
            get_samplesheet_records(samplesheet_list)
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
        merged_sample_data.fillna(0, inplace=True)
        (labels, total_cluster_raw, total_cluster_pf, total_yield) = \
            get_flowcell_summary_data(summary_data=sum_df)
        flow_cell_data = {
            'labels': labels,
            'total_cluster_raw': total_cluster_raw,
            'total_cluster_pf': total_cluster_pf,
            'total_yield': total_yield}
        lane_plots = \
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
            get_flowcell_project_summary_plot_for_db(
                summary_data=sum_df,
                sample_data=merged_sample_data)
        undetermined_plots = \
            get_undetermined_plot(
                undetermined_data=undetermined_data)
        undetermined_tables = \
            get_undetermined_table(
                undetermined_data=undetermined_data)
        json_data = {
            'run_name': run_name,
            'samplesheet_tag': samplesheet_tag,
            'flowcell_cluster_plot': flow_cell_data,
            'project_summary_table': project_summary_data,
            'project_summary_plot': project_summary_plot,
            'sample_table': table_data,
            'sample_plot': lane_plots,
            'undetermined_table': undetermined_tables,
            'undetermined_plot': undetermined_plots}
        with tempfile.TemporaryDirectory() as temp_dir :
            temp_json_file = \
                os.path.join(
                    temp_dir,
                    '{0}_{1}.json'.format(run_name, samplesheet_tag))
            with open(temp_json_file, 'w') as fp:
                json.dump(json_data, fp)
            copy2(temp_json_file, final_json_file)
    except Exception as e:
        raise ValueError("Failed to create json output, error: {0}".format(e))
