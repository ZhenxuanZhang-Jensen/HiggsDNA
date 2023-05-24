import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
eff=pd.read_csv("/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/SLm1200/UL17_R_gghh_SL_M-1200_2017/cutflow_eff.csv")
df=pd.read_csv("/eos/user/s/shsong/combined_WWgg/parquet/sig_latest/SLm1200/UL17_R_gghh_SL_M-1200_2017/event_yield.csv")

def create_event_selection_table(numdataframe,effdataframe):
    new_columns = {}
    for column in numdataframe.columns:
        if 'event number' in column:
            new_column = column.replace('event number', '')
            new_columns[column] = new_column
    
    numdataframe = numdataframe.rename(columns=new_columns)
    
    # 提取变量名包含"event"的列
    event_columns = [col for col in numdataframe.columns if 'event' in col]
    
    # 创建包含这些列的新数据框
    event_data = numdataframe[event_columns]
    event_data_transposed = event_data.transpose()
    event_columns = [col for col in effdataframe.columns if 'event efficiency' in col]
    event_eff=effdataframe[event_columns]
    new_columns = {}
    for column in event_eff.columns:
        if 'event efficiency' in column:
            new_column = column.replace('event efficiency', '')
            new_columns[column] = new_column
        
    event_eff = event_eff.rename(columns=new_columns)
    # # event_eff.transpose()
    eff_columns = [col for col in event_eff.columns if 'event:' in col]
    eventeff=event_eff[eff_columns]
    eventeff_transposed=eventeff.transpose()
    return event_data_transposed,eventeff_transposed
eventnum,eventeff=create_event_selection_table(df,eff)
