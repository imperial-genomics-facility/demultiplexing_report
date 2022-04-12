import os
import subprocess
import numpy as np
import pandas as pd
from io import StringIO
from iplotter import ChartJSPlotter
from iplotter import GCPlotter

def combine_bclconvert_demultiplex_stats_csv(demultiplex_stats_list: list) \
        -> pd.DataFrame:
    try:
        merged_df = pd.DataFrame()
        for entry in demultiplex_stats_list:
            if not os.path.exists(entry):
                raise IOError('Missing file {0}'.format(entry))
            df = pd.read_csv(entry)
            if len(merged_df.index) == 0:
                merged_df = df.copy()
            else:
                merged_df = \
                    pd.concat(
                        [merged_df, df],
                        ignore_index=True)
        expected_columns = [
            'Lane', 'SampleID', 'Index', '# Reads',
            '# Perfect Index Reads', '# One Mismatch Index Reads',
            '# Two Mismatch Index Reads', '% Reads',
            '% Perfect Index Reads', '% One Mismatch Index Reads',
            '% Two Mismatch Index Reads']
        for i in expected_columns:
            if i not in merged_df.columns:
                raise KeyError(
                    "Missing column {0} in demultiplex_stats".\
                    format(i))
        combined_df = \
            merged_df.\
            groupby([
                'Lane', 'SampleID', 'Index']).\
            agg({
                '# Reads': np.sum,
                '# Perfect Index Reads': np.sum,
                '# One Mismatch Index Reads': np.sum,
                '# Two Mismatch Index Reads': np.sum,
                '% Reads': np.mean,
                '% Perfect Index Reads': np.mean,
                '% One Mismatch Index Reads': np.mean,
                '% Two Mismatch Index Reads': np.mean})
        combined_df.reset_index(inplace=True)
        return combined_df
    except Exception as e:
        raise ValueError(
                'Failed to combine Demultiplex_Stats csv files, error: {0}'.\
                format(e))

def combine_bclconvert_quality_metrics_csv(quality_metrics_list: list) \
        -> pd.DataFrame:
    try:
        merged_df = pd.DataFrame()
        for entry in quality_metrics_list:
            if not os.path.exists(entry):
                raise IOError('Missing file {0}'.format(entry))
            df = pd.read_csv(entry)
            if len(merged_df.index) == 0:
                merged_df = df.copy()
            else:
                merged_df = \
                    pd.concat(
                        [merged_df, df],
                        ignore_index=True)
        expected_columns = [
            "Lane", "SampleID", "index", "index2", "ReadNumber",
            "Yield", "YieldQ30", "QualityScoreSum", "% Q30"]
        for i in expected_columns:
            if i not in merged_df.columns:
                raise KeyError(
                    "Missing column {0} in quality_metrics".\
                    format(i))
        combined_df = \
            merged_df.\
            groupby([
                "Lane", "SampleID", "index", "index2", "ReadNumber"]).\
            agg({
                "Yield": np.sum,
                "YieldQ30": np.sum,
                "QualityScoreSum": np.sum,
                "Mean Quality Score (PF)": np.mean,
                "% Q30": np.mean})
        combined_df.reset_index(inplace=True)
        return combined_df
    except Exception as e:
        raise ValueError(
                'Failed to combine Quality_Metrics csv files, error: {0}'.\
                format(e))

def combine_bclconvert_top_unknown_barcodes_csv(
        top_unknown_barcodes_list: list) \
        -> pd.DataFrame:
    try:
        merged_df = pd.DataFrame()
        for entry in top_unknown_barcodes_list:
            if not os.path.exists(entry):
                raise IOError('Missing file {0}'.format(entry))
            df = pd.read_csv(entry)
            if len(merged_df.index) == 0:
                merged_df = df.copy()
            else:
                merged_df = \
                    pd.concat(
                        [merged_df, df],
                        ignore_index=True)
        expected_columns = [
            'Lane', 'index', 'index2', '# Reads',
            '% of Unknown Barcodes', '% of All Reads']
        for i in expected_columns:
            if i not in merged_df.columns:
                raise KeyError(
                    "Missing column {0} in quality_metrics".\
                    format(i))
        combined_df = \
            merged_df.\
            groupby([
                'Lane', 'index', 'index2']).\
            agg({
                '# Reads': np.sum,
                '% of Unknown Barcodes': np.mean,
                '% of All Reads': np.mean})
        combined_df.reset_index(inplace=True)
        combined_df = \
            combined_df.\
            sort_values('# Reads', ascending=False).\
            head(20).\
            sort_values([
                'Lane', '# Reads'], 
                ascending=False)
        return combined_df
    except Exception as e:
        raise ValueError(
                'Failed to combine Top_Unknown_Barcodes csv files, error: {0}'.\
                format(e))

def combine_bclconvert_index_hopping_counts_csv(
        index_hopping_counts_list: list) \
        -> pd.DataFrame:
    try:
        merged_df = pd.DataFrame()
        for entry in index_hopping_counts_list:
            if not os.path.exists(entry):
                raise IOError('Missing file {0}'.format(entry))
            df = pd.read_csv(entry)
            if len(merged_df.index) == 0:
                merged_df = df.copy()
            else:
                merged_df = \
                    pd.concat(
                        [merged_df, df],
                        ignore_index=True)
        expected_columns = [
            'Lane', 'SampleID', 'index', 'index2', '# Reads',
            '% of Hopped Reads', '% of All Reads']
        for i in expected_columns:
            if i not in merged_df.columns:
                raise KeyError(
                    "Missing column {0} in quality_metrics".\
                    format(i))
        merged_df['SampleID'] = \
            merged_df['SampleID'].fillna('UNKNOWN')
        merged_df.dropna(inplace=True)
        filt_merged_df = \
            merged_df[merged_df['% of Hopped Reads'] > 0.05]
        combined_df = \
            filt_merged_df.\
            groupby([
                'Lane', 'SampleID', 'index', 'index2']).\
            agg({
                '# Reads': np.sum,
                '% of Hopped Reads': np.mean,
                '% of All Reads': np.mean})
        combined_df.reset_index(inplace=True)
        return combined_df
    except Exception as e:
        raise ValueError(
                'Failed to combine Index_Hopping_Counts csv files, error: {0}'.\
                format(e))

def get_samplesheet_records(samplesheets: list) \
        -> pd.DataFrame:
    try:
        all_samplesheet_data = pd.DataFrame()
        for f in samplesheets:
            data_section = False
            samplesheet_data_list = list()
            with open(f, 'r') as fp:
                for i in fp:
                    i = i.strip()
                    if i == '':
                        continue
                    if i.startswith('['):
                        data_section = False
                        if i.startswith('[Data]') or \
                           i.startswith('[data]') or \
                           i.startswith('[BCLConvert_Data]'):
                            data_section = True
                            continue
                    if data_section:
                        samplesheet_data_list.\
                            append(i.split(','))
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
        raise ValueError(e)

def get_interop_index_stats(run_dir: str) \
        -> pd.DataFrame:
    try:
        index_csv = \
            subprocess.\
            check_output([
                "interop_index-summary",
                run_dir,
                "--csv=1"])
        if isinstance(index_csv, bytes):
            index_csv = index_csv.decode('utf-8')
        index_csv = index_csv.split("\n")
        counter = 99
        key = None
        lane_data = list()
        total_reads = dict()
        for i in index_csv:
            if i.startswith('Lane'):
                counter = 0
                key = i.split(",")[0]
                if key is not None:
                    key = key.replace("Lane", "").strip()
                lane_data = list()
                continue
            if counter < 2:
                lane_data.append(i)
                counter += 1
            if counter == 2:
                total_reads.update({key: lane_data})
        formatted_lane_data = dict()
        for lane_id, lane_data in total_reads.items():
            csv_data = StringIO('\n'.join(total_reads.get(lane_id)))
            data = pd.read_csv(csv_data).to_dict(orient='records')
            formatted_lane_data.update({lane_id: data})
        flowcell_summary_data = list()
        for lane_id, lane_data in formatted_lane_data.items():
            flowcell_summary_data.\
                append({
                    'Lane': lane_id,
                    'Total Reads': lane_data[0].get('Total Reads'),
                    'PF Reads': lane_data[0].get('PF Reads')})
        flowcell_summary_data = \
            pd.DataFrame(flowcell_summary_data)
        return flowcell_summary_data
    except Exception as e:
        raise ValueError(e)

def merge_known_samples(
    demultiplex_stats_df: pd.DataFrame,
    quality_metrics_df: pd.DataFrame,
    samplesheet_df: pd.DataFrame) \
        -> pd.DataFrame:
    try:
        temp_demultiplex_stats_df = \
            demultiplex_stats_df.\
            copy()
        temp_quality_metrics_df = \
            quality_metrics_df.\
            copy()
        temp_samplesheet_df = \
            samplesheet_df.\
            copy()
        filt_quality_metrics_df = \
            temp_quality_metrics_df[
                ~temp_quality_metrics_df['ReadNumber'].
                str.startswith("I")]
        if 'index2' in filt_quality_metrics_df.columns:
            filt_quality_metrics_df['Index'] = \
                filt_quality_metrics_df[['index', 'index2']].\
                agg('-'.join, axis=1)
        else:
            filt_quality_metrics_df['Index'] = \
                filt_quality_metrics_df['index'].copy()
        agg_combined_qmetrics_df = \
            filt_quality_metrics_df.\
            groupby(['Lane', 'SampleID', 'Index']).\
            agg({
                'Yield': 'sum',
                'Mean Quality Score (PF)': 'mean',
                '% Q30': 'mean'})
        joined_data1 = \
            temp_demultiplex_stats_df.\
            set_index(['Lane', 'SampleID', 'Index']).\
            join(agg_combined_qmetrics_df, how='left').\
            reset_index()
        joined_data1 = \
            joined_data1[[
                'Lane',
                'SampleID',
                'Index',
                '# Reads',
                '% Reads',
                '% Perfect Index Reads',
                'Yield',
                '% Q30',
                'Mean Quality Score (PF)']]
        joined_data1.\
            set_index('SampleID', inplace=True)
        temp_samplesheet_df.\
            set_index('Sample_ID', inplace=True)
        final_joined_data = \
            joined_data1.\
            join(temp_samplesheet_df[['Sample_Name', 'Sample_Project']].\
                 drop_duplicates(),
                 how='left')
        final_joined_data.index.rename('SampleID', inplace=True)
        final_joined_data.reset_index(inplace=True)
        final_joined_data = \
            final_joined_data[[
                'Lane',
                'Sample_Project',
                'SampleID',
                'Sample_Name',
                'Index',
                '# Reads',
                '% Reads',
                '% Perfect Index Reads',
                'Yield',
                '% Q30',
                'Mean Quality Score (PF)']]
        final_joined_data['Sample_Project'].\
            fillna('UNKNOWN', inplace=True)
        final_joined_data['Sample_Name'].\
            fillna('UNKNOWN', inplace=True)
        final_joined_data.\
            sort_values([
                'Lane',
                'Sample_Project',
                '# Reads'],
                ascending=[True, True, False],
                inplace=True)
        return final_joined_data
    except Exception as e:
        raise ValueError(e)

def get_flowcell_summary_plot(flowcell_summary_df: pd.DataFrame) -> dict:
    try:
        temp_flowcell_summary_df = \
            flowcell_summary_df.\
            copy()
        labels = \
            temp_flowcell_summary_df['Lane'].\
            map(lambda x: 'Lane {0}'.format(x)).\
            values.tolist()
        datasets = [{
            "label": "Total cluster raw",
            "data": temp_flowcell_summary_df["Total Reads"].astype(int).values.tolist(),
            "backgroundColor": 'rgba(255, 99, 132, 0.8)',
            "borderColor": 'rgba(255, 99, 132, 0.8)',
            "borderWidth": 1},{
            "label": "Total cluster pf",
            "data": temp_flowcell_summary_df["PF Reads"].astype(int).values.tolist(),
            "backgroundColor": 'rgba(54, 162, 235, 0.8)',
            "borderColor": 'rgba(54, 162, 235, 0.8)',
            "borderWidth": 1}]
        options = {
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
        data = {
            "labels": labels,
            "datasets": datasets
        }
        plot = \
            cj_plotter.\
            plot(
                data,
                'bar',
                options=options,
                w=400,
                h=300)
        return plot
    except Exception as e:
        raise ValueError(e)

def get_flowcell_project_summary_plot(merged_data: pd.DataFrame) \
        -> dict:
    try:
        temp_merged_data = \
            merged_data.\
            copy()
        project_groups = \
            temp_merged_data.\
            groupby(['Lane', 'Sample_Project']).\
            agg({'# Reads': 'sum'}).\
            reset_index()
        project_group_df = list()
        for lane_id, l_data in project_groups.groupby('Lane'):
            data_row = {'Lane': 'Lane {0}'.format(lane_id)}
            for project_id, p_data in l_data.groupby('Sample_Project'):
                data_row.update({project_id: p_data['# Reads'].sum()})
            project_group_df.append(data_row)
        project_group_df = \
            pd.DataFrame(project_group_df)
        data = list()
        data.append(project_group_df.columns.tolist())
        data.extend(project_group_df.values.tolist())
        options = {
            "title": 'Project summary plot',
            "width": 600,
            "height": 400,
            "chartArea": {"left": 50, "width": "60%"},
            "legend": {"position": 'right', "maxLines": 20, "fontSize": 3},
            "dataOpacity": 0.5,
            "colors": [
                '#0173B2', '#DE8F05', '#029E73', '#D55E00', '#CC78BC', '#CA9161',
                '#FBAFE4', '#949494', '#ECE133', '#56B4E9', '#0173B2'],
            "bar": {"groupWidth": '70%'},
            "isStacked": "percent"}
        gcplotter = GCPlotter()
        plot = \
            gcplotter.\
            plot(
                data,
                chart_type="BarChart",
                chart_package='corechart',
                options=options)
        return plot
    except Exception as e:
        raise ValueError(e)

def get_flowcell_project_summary_table(merged_data: pd.DataFrame) \
        -> pd.DataFrame:
    try:
        temp_merged_data = \
            merged_data.\
            copy()
        if 'Lane' in temp_merged_data.columns:
            summary_df = \
                temp_merged_data.\
                groupby(['Lane', 'Sample_Project']).\
                agg(len)['SampleID'].\
                reset_index()
        else:
            summary_df = \
                temp_merged_data.\
                groupby(['Sample_Project',]).\
                agg(len)['SampleID'].\
                reset_index()
        return summary_df
    except Exception as e:
        raise ValueError(e)

def get_sample_dist_plots(
        merged_data: pd.DataFrame,
        get_plots: bool = True) \
        -> dict:
    try:
        temp_merged_data = \
            merged_data.\
            copy()
        lane_plots = dict()
        cj_plotter = ChartJSPlotter()
        bg_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        border_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        options = {
            'scales': {
                'y': {
                    'beginAtZero': True
                },
                'xAxes': [{
                    'ticks': {
                        'fontSize': 8
                    }
                }]
            },
            'responsive': True,
            'plugins': {
                'legend': {
                    'position': 'top'}
            }}
        if 'Lane' in temp_merged_data.columns:
            for lane_id, l_data in temp_merged_data.groupby('Lane'):
                lane_samples = l_data['SampleID'].values.tolist()
                datasets = list()
                counter = 0
                for project_name, p_data in l_data.groupby('Sample_Project'):
                    read_counts = \
                        p_data.\
                        set_index('SampleID')['# Reads'].\
                        reindex(lane_samples).\
                        fillna(0).\
                        values.\
                        tolist()
                    datasets.append({
                        "label": project_name,
                        "data": read_counts,
                        "backgroundColor": bg_colors[counter],
                        "borderColor": border_colors[counter]})
                    counter += 1
                data = {
                    "labels": lane_samples,
                    "datasets": datasets}
                if get_plots:
                    plot = \
                        cj_plotter.\
                        plot(
                            data,
                            'bar',
                            options=options,
                            w=800,
                            h=400)
                    lane_plots.update({int(lane_id): plot})
                else:
                    lane_plots.update({int(lane_id): data})
        else:
            lane_samples = temp_merged_data['SampleID'].values.tolist()
            datasets = list()
            counter = 0
            for project_name, p_data in temp_merged_data.groupby('Sample_Project'):
                read_counts = \
                    p_data.\
                    set_index('SampleID')['# Reads'].\
                    reindex(lane_samples).\
                    fillna(0).\
                    values.\
                    tolist()
                datasets.append({
                    "label": project_name,
                    "data": read_counts,
                    "backgroundColor": bg_colors[counter],
                    "borderColor": border_colors[counter]})
                counter += 1
            data = {
                "labels": lane_samples,
                "datasets": datasets}
            if get_plots:
                plot = \
                    cj_plotter.\
                    plot(
                        data,
                        'bar',
                        options=options,
                        w=800,
                        h=400)
                lane_plots.update({1: plot})
            else:
                lane_plots.update({1: data})
        return lane_plots
    except Exception as e:
        raise ValueError(e)

def get_undetermined_table(
        unknown_df: pd.DataFrame) \
        -> pd.DataFrame:
    try:
        temp_unknown_df = \
            unknown_df.\
            copy()
        if 'index2' in temp_unknown_df.columns:
            temp_unknown_df['Barcode'] = \
                temp_unknown_df[['index', 'index2']].\
                agg('-'.join, axis=1)
            temp_unknown_df['Barcode_I1_RC'] = \
                temp_unknown_df['index'].\
                str.translate(str.maketrans('ATGC', 'TACG'))[::-1] + \
                '-' + temp_unknown_df['index2']
            temp_unknown_df['Barcode_I2_RC'] = \
                temp_unknown_df['index'] + '-' + \
                temp_unknown_df['index2'].\
                str.translate(str.maketrans('ATGC', 'TACG'))[::-1]
        else:
            temp_unknown_df['Barcode'] = \
                temp_unknown_df['index'].\
                copy()
            temp_unknown_df['Barcode_I1_RC'] = \
                temp_unknown_df['index'].\
                str.translate(str.maketrans('ATGC', 'TACG'))[::-1]
            temp_unknown_df['Barcode_I2_RC'] = ''
        temp_unknown_df = \
            temp_unknown_df[[
                'Lane',
                '# Reads',
                'Barcode',
                'Barcode_I1_RC',
                'Barcode_I2_RC']]
        temp_unknown_df.\
            sort_values([
                'Lane',
                '# Reads'],
                ascending=[True, False],
                inplace=True)
        return temp_unknown_df
    except Exception as e:
        raise ValueError(e)


def get_undetermined_plots(
        unknown_df: pd.DataFrame,
        get_plots: bool = True) \
        -> dict:
    try:
        bg_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        border_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        cj_plotter = ChartJSPlotter()
        temp_unknown_df = \
            unknown_df.\
            copy()
        if 'index2' in temp_unknown_df.columns:
            temp_unknown_df['Index'] = \
                temp_unknown_df[['index', 'index2']].\
                agg('-'.join, axis=1)
        else:
            temp_unknown_df['Index'] = \
                temp_unknown_df['index'].\
                copy()
        undetermined_plots = dict()
        options = {
            'scales': {
                'y': {
                    'beginAtZero': True
                },
                'xAxes': [{
                    'ticks': {
                        'fontSize': 8
                    }
                }]
            },
            'responsive': True,
            'plugins': {
                'legend': {
                    'position': 'top'}
            }}
        if 'Lane' in temp_unknown_df.columns:
            for lane_id, l_data in temp_unknown_df.groupby('Lane'):
                barcode_labels = l_data['Index'].values.tolist()
                barcode_count = l_data['# Reads'].values.tolist()
                data = {
                    "labels": barcode_labels,
                    "datasets": [{
                        'label': 'Undetermined - Lane {0}'.format(lane_id),
                        'data': barcode_count,
                        "backgroundColor": bg_colors[0],
                        "borderColor": border_colors[0]}]}
                if get_plots:
                    plot = \
                        cj_plotter.\
                        plot(
                            data,
                            'bar',
                            options=options,
                            w=800,
                            h=400)
                    undetermined_plots.\
                        update({int(lane_id): plot})
                else:
                    undetermined_plots.\
                        update({int(lane_id): data})
        else:
            barcode_labels = temp_unknown_df['Index'].values.tolist()
            barcode_count = temp_unknown_df['# Reads'].values.tolist()
            data = {
                "labels": barcode_labels,
                "datasets": [{
                    'label': 'Undetermined - Lane {0}'.format(1),
                    'data': barcode_count,
                    "backgroundColor": bg_colors[0],
                    "borderColor": border_colors[0]}]}
            if get_plots:
                plot = \
                    cj_plotter.\
                    plot(
                        data,
                        'bar',
                        options=options,
                        w=800,
                        h=400)
                undetermined_plots.\
                    update({1: plot})
            else:
                undetermined_plots.\
                    update({1: data})
        return undetermined_plots
    except Exception as e:
        raise ValueError(e)


def get_hop_plot(
        hop_df: pd.DataFrame,
        get_plot: bool = True) \
        -> dict:
    try:
        lane_data = list()
        options = {
            'scales': {
                'y': {
                    'beginAtZero': True
                },
                'xAxes': [{
                    'ticks': {
                        'fontSize': 8
                    }
                }]
            },
            'responsive': True,
            'plugins': {
                'legend': {
                    'position': 'top'}
            }}
        bg_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        border_colors = [
            'rgba(255, 99, 132, 0.8)',
            'rgba(54, 162, 235, 0.8)',
            'rgba(255, 206, 86, 0.8)',
            'rgba(75, 192, 192, 0.8)',
            'rgba(153, 102, 255, 0.8)',
            'rgba(255, 159, 64, 0.8)',
            'rgba(255, 159, 10, 0.8)',
            'rgba(255, 159, 192, 0.8)']
        cj_plotter = ChartJSPlotter()
        temp_hop_df = hop_df.copy()
        if 'index2' in temp_hop_df:
            temp_hop_df['Index'] = \
                temp_hop_df[['index', 'index2']].\
                agg('-'.join, axis=1)
        else:
            temp_hop_df['Index'] = \
                temp_hop_df['index'].copy()
        if 'Lane' in temp_hop_df:
            for lane_id, l_data in temp_hop_df.groupby('Lane'):
                data_row = {'Lane': lane_id}
                for index, i_data in l_data.groupby('Index'):
                    data_row.update({index: i_data['# Reads'].sum()})
                lane_data.append(data_row)
        else:
            data_row = {'Lane': 1}
            for index, i_data in temp_hop_df.groupby('Index'):
                data_row.update({index: i_data['# Reads'].sum()})
            lane_data.append(data_row)
        lane_data = pd.DataFrame(lane_data)
        lane_data.fillna(0, inplace=True)
        barcodes = [c for c in lane_data.columns if c != 'Lane']
        dataset = list()
        counter = 0
        for lane_id, l_data in lane_data.groupby('Lane'):
            barcode_count = l_data[barcodes].values.tolist()[0]
            dataset.append({
                'label': 'Hopping - Lane {0}'.format(lane_id),
                'data': barcode_count,
                "backgroundColor": bg_colors[counter],
                "borderColor": border_colors[counter]})
            counter += 1
        data = {
            'labels': barcodes,
            'datasets': dataset}
        if get_plot:
            data = \
                cj_plotter.\
                plot(
                    data,
                    'bar',
                    options=options,
                    w=800,
                    h=400)
        return data
    except Exception as e:
        raise ValueError(e)

def get_demult_report_and_plots_for_bclconvert(
        run_dir: str,
        reports_dir: str) -> None:
    try:
        if not os.path.exists(run_dir):
            raise IOError('Missing run dir {0}'.format(run_dir))
        required_report_csvs = [
            'Demultiplex_Stats.csv',
            'Quality_Metrics.csv',
            'Top_Unknown_Barcodes.csv',
            'SampleSheet.csv']
        for f in required_report_csvs:
            filepath = \
                os.path.join(reports_dir, f)
            if not os.path.exists(filepath):
                raise IOError('Missing report file {0}'.format(filepath))
        combined_demux_df = \
            combine_bclconvert_demultiplex_stats_csv([
                os.path.join(reports_dir, 'Demultiplex_Stats.csv')])
        combined_qmetrics_df = \
            combine_bclconvert_quality_metrics_csv([
                os.path.join(reports_dir, 'Quality_Metrics.csv')])
        combined_unknown_df = \
            combine_bclconvert_top_unknown_barcodes_csv([
                os.path.join(reports_dir, 'Top_Unknown_Barcodes.csv')])
        combined_ihop_df = pd.DataFrame()
        hop_plot = {}
        if os.path.exists(os.path.join(reports_dir, 'Index_Hopping_Counts.csv')):
            combined_ihop_df = \
                combine_bclconvert_index_hopping_counts_csv([
                    os.path.join(reports_dir, 'Index_Hopping_Counts.csv')])
            if len(combined_ihop_df.index) > 0:
                hop_plot = \
                    get_hop_plot(combined_ihop_df)
        samplesheet_df = \
            get_samplesheet_records(
                samplesheets=[os.path.join(reports_dir, 'SampleSheet.csv')])
        flowcell_summary_data = \
            get_interop_index_stats(run_dir=run_dir)
        flowcell_summary_data_plot = \
            get_flowcell_summary_plot(flowcell_summary_data)
        merged_df = \
            merge_known_samples(
                demultiplex_stats_df=combined_demux_df,
                quality_metrics_df=combined_qmetrics_df,
                samplesheet_df=samplesheet_df)
        flowcell_project_summary_plot = \
            get_flowcell_project_summary_plot(merged_df)
        flowcell_project_summary_table = \
            get_flowcell_project_summary_table(merged_df)
        sample_dist_plots = \
            get_sample_dist_plots(merged_df)
        undetermined_plots = \
            get_undetermined_plots(combined_unknown_df)
        undetermined_table = \
            get_undetermined_table(combined_unknown_df)
        return flowcell_summary_data_plot, flowcell_project_summary_plot, \
            merged_df, flowcell_project_summary_table, sample_dist_plots, \
            undetermined_plots, undetermined_table, combined_ihop_df, \
            hop_plot
    except Exception as e:
        raise ValueError(e)