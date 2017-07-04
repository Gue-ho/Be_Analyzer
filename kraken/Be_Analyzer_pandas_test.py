import os
import pytest
import pandas

####### fuction

def Seq_standard(wt_seq, target_seq, filt_r,end_range):
    for i in range(22,len(wt_seq)-23):
        if wt_seq[i:i+23] == target_seq:
            start_pos=i+17-end_range
            end_pos=i+17+end_range
            
            if start_pos < 0:
                start_pos=0
            if end_pos > len(wt_seq):
                end_pos = len(wt_seq)
                
            s_seq = wt_seq[i+17-filt_r:i+17+filt_r]
            seq_range = wt_seq[start_pos:end_pos]
            break
        
        elif wt_seq[i:i+23] == target_seq.translate(target_seq.maketrans('ATGC','TACG'))[::-1]:
            start_pos=i+6-end_range
            end_pos=i+6+end_range
            
            if start_pos < 0:
                start_pos=0
            if end_pos > len(wt_seq):
                end_pos = len(wt_seq)
                
            s_seq = wt_seq[i+6-filt_r:i+6+filt_r]
            seq_range = wt_seq[start_pos:end_pos]
            break
        else:
            start_pos=0
            end_pos=len(wt_seq)
            
            s_seq = 'None'
            seq_range = wt_seq

    return seq_range, s_seq

def Indicator_table(seq_range, indicator_range):
    pri_for = seq_range.upper()[:indicator_range]
    pri_back = seq_range.upper()[-indicator_range:]                 
    length_range=len(seq_range)

    for_table=[pri_for]
    for i in range(len(pri_for)):
        for n in ['A','T','G','C']:
            for_table.append(pri_for[:i]+n+pri_for[i+1:])

    back_table=[pri_back]
    for i in range(len(pri_back)):
        for n in ['A','T','G','C']:
            back_table.append(pri_back[:i]+n+pri_back[i+1:])

    return for_table, back_table

def Seq_check(out_seq,for_table,back_table):
    for_check=0
    back_check=0

    if out_seq.find('N')!=-1:
        return 0
    
    for pri in for_table:
        if out_seq.find(pri)!=-1:
            for_check=1
            start_seq=out_seq.find(pri)

    for back in back_table:
        if out_seq.find(back)!=-1:
            back_check=1
            end_seq=out_seq.find(back)+len(back_table[0])

    if for_check==1 and back_check==1:
        acu_seq=out_seq[start_seq:end_seq]
        return acu_seq
    else:
        return 0

def Read_seq(file_name,for_table,back_table,filt_n):#Seq_check

    all_count=0
    primer_count=0
    filted_count=0

    line_list={}
    
    f=open(file_name).readlines()
    
    for n in range(int(len(f)/4)):

        out_seq=f[4*n+1]
        all_count+=1

        acu_seq=Seq_check(out_seq,for_table,back_table)
        
        if acu_seq==0:
            continue

        primer_count+=1
        
        try:
            line_list[acu_seq]+=1
            continue
        except:
            line_list[acu_seq]=1

    del_dict=[]
    
    for x,y in line_list.items():
        if y<=filt_n:
            filted_count+=1
            del_dict.append(x)

    for x in del_dict:
        del line_list[x]

    return line_list, all_count, primer_count, filted_count

def Substitution(emboss_wt,emboss_seq):
    sub_wt=''
    sub_seq=''
    sub_sym=''
    
    for i in range(len(emboss_wt)):
        if emboss_wt[i]!='-' and emboss_seq[i]!='-':
            if emboss_wt[i]!=emboss_seq[i]:
                sub_wt+=emboss_wt[i]
                sub_seq+=emboss_seq[i].lower()
                sub_sym+='.'
            else:
                sub_wt+=emboss_wt[i]
                sub_seq+=emboss_seq[i]
                sub_sym+='|'
        else:
            sub_wt+=emboss_wt[i]
            sub_seq+=emboss_seq[i]
            sub_sym+=' '

    return sub_wt, sub_seq, sub_sym

def Indel_count(emboss_wt,emboss_seq,sub_wt,sub_seq,count,Count_list,Substitution_position,Substitution_pattern,CtoD):
    if emboss_wt.find('-')!=-1 and emboss_seq.find('-')!=-1:
        Indel='Others'
        Count_list[5]+=count
    elif emboss_wt.find('-')!=-1:
        Indel='Insertion'
        Count_list[3]+=count
    elif emboss_seq.find('-')!=-1:
        Indel='Deletion'
        Count_list[4]+=count
    else:
        if sub_seq.upper()!=sub_seq:
            Indel='Substituion'
            Count_list[1]+=count
            for i in range(len(sub_seq)):
                if sub_seq[i] in ['a','t','g','c']:
                    Substitution_position[i]+=count
                    Substitution_pattern[i]['atgc'.find(sub_seq[i])]+=count
                    if sub_wt[i]=='C':
                        Count_list[2]+=count
                        CtoD[i]['atgc'.find(sub_seq[i])]+=count
        else:
            Indel='WT'
            Count_list[0]+=count
    return Count_list, Indel, Substitution_position, Substitution_pattern, CtoD
    

def Read_EMBOSS(file, line_list, seq_range):#Substitution, Indel_count

    Count_list=[0,0,0,0,0,0] # [WT_count,Substituion_count,CtoD,Insertion_count,Deletion_count,Others]
    emboss_dict={}

    Substitution_position=[]
    Substitution_pattern=[]
    CtoD=[]

    for i in range(len(seq_range)): #[A,T,G,C]
        Substitution_position.append(0)
        Substitution_pattern.append([0,0,0,0])
        CtoD.append([0,0,0,0])
    
    f2=open(file).readlines()

    emboss_wt=''
    emboss_seq=''
    for i in range(len(f2)):
        if len(f2[i].split())>1:
            if f2[i].split()[0]=='EMBOSS_001':
                emboss_wt+=f2[i].split()[2]
                emboss_seq+=f2[i+2].split()[2]

            if f2[i].split()[0]=='EMBOSS_001' and f2[i+4]=='\n':
                ori_seq=emboss_seq.replace('-','')
                
                (sub_wt,sub_seq,sub_sym) = Substitution(emboss_wt,emboss_seq)

                (Count_list, Indel, Substitution_position, Substitution_pattern, CtoD)=Indel_count(emboss_wt,emboss_seq,sub_wt,sub_seq,line_list[ori_seq],Count_list,Substitution_position,Substitution_pattern,CtoD)
                
                emboss_dict[ori_seq]=[sub_wt,sub_sym,sub_seq,line_list[ori_seq],Indel]

                emboss_wt=''
                emboss_seq=''
        

    return emboss_dict, Count_list, Substitution_position, Substitution_pattern, CtoD

def writer(line_list):
    fw=open('result.txt','w')
    for x, y in line_list.items():
        fw.write(emboss_list[x][0]+'\t'+str(y)+'\n')
        fw.write(emboss_list[x][1]+'\n')
        fw.write(emboss_list[x][2]+'\n\n')
    fw.close()
    return 0

def emboss_needle(seq_range, line_list, direc):

    if direc[-1]!='/':
        direc+='/'
        
    f1=open('wt_seq_emboss.txt','w')
    f1.write('>\n'+seq_range)
    f1.close()

    f2=open('emboss_before.txt','w')
    for x,y in line_list.items():
        f2.write('>\n'+x+'\n')
    f2.close()

    os.system('needle -asequence '+direc+'/wt_seq_emboss.txt -bsequence '+direc+'/emboss_before.txt -gapopen 10 -gapextend 0.5 -outfile '+direc+'/emboss_after.txt')

    return 0

def Total_Substitutions(seq_range,target_seq,Substitution_position,Substitution_pattern):
    wb1=pandas.ExcelWriter('Total_Substitutions.xlsx',engine='xlsxwriter')
    df=pandas.DataFrame(columns=('WT_seq','Substitution_count','A','T','G','C'))
    row=0
    for i in range((seq_range.find(target_seq)+17)*(-1),(seq_range.find(target_seq)+17)*(-1)+len(seq_range)+1):
        if i==0:
            continue
        row+=1
        df.loc[i]=[seq_range[row-1],Substitution_position[row-1]]+Substitution_pattern[row-1]
    df.columns.name='Site'
    df.to_excel(wb1,'Sheet1')
    wb1.save()
    wb1.close()
    
    return 0

def C_to_D_Substitutions(seq_range,target_seq,CtoD):
    wb2=pandas.ExcelWriter('C_to_D_Substitutions.xlsx',engine='xlsxwriter')
    df=pandas.DataFrame(columns=('WT_seq','substitution','A','T','G','C'))
    df.columns.name='Site'
    row=0
    for i in range((seq_range.find(target_seq)+17)*(-1),(seq_range.find(target_seq)+17)*(-1)+len(seq_range)+1):
        if i==0:
            continue
        row+=1
        all_c=0
        for x in CtoD[row-1]:
            all_c+=x
        df.loc[i]=[seq_range[row-1],all_c]+CtoD[row-1]
    df.to_excel(wb2,'Sheet')
    wb2.save()
    wb2.close()

    return 0

def Transition_figure(seq_range,Substitution_count,Substitution_pattern):
    mat=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    for i in range(len(seq_range)):
        x=int(seq_range[i].translate(seq_range[i].maketrans('ATGC','0123')))
        for y in range(4):
            mat[x][y]+=Substitution_pattern[i][y]
    a_mat=[0,0,0,0]
    for i in range(len(mat)):
        for x in mat[i]:
            a_mat[i]+=x
    wb3=pandas.ExcelWriter('Substitution_transition.xlsx',engine='xlsxwriter')
    df=pandas.DataFrame(columns=('A','T','G','C'))
    res_mat=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    for x in range(4):
        for y in range(4):
            if a_mat[x]!=0:
                res_mat[y][x]=str(round(mat[x][y]*100/a_mat[x],5))
            elif a_mat[x]==0:
                res_mat[y][x]=str(0)
    for i in ['A','T','G','C']:
        x=int(i.translate(i.maketrans('ATGC','0123')))
        df.loc[i]=res_mat[x]
    df.to_excel(wb3,sheet_name='Sheet1',startcol=1,startrow=1)
    workbook=wb3.book
    ws3=wb3.sheets['Sheet1']
    cell_format=workbook.add_format({'align':'center','valign':'vcenter','bold':True})
    ws3.merge_range('A3:A6','',cell_format)
    ws3.merge_range('C1:F1','',cell_format)
    ws3.write_rich_string('A3','Mutation',cell_format)
    ws3.write_rich_string('C1','WT',cell_format)
    ws3.write_row(1,2,['A','T','G','C'],cell_format)
    ws3.write_column(2,1,['A','T','G','C'],cell_format)
    wb3.save()
    wb3.close()

    return 0

def Highest(vals):
    res=vals[0]
    c=0
    for i in range(len(vals)):
        if res<vals[i]:
            res=vals[i]
            c=i
    return c

def Chart(seq_range,target_seq,all_c,Substitution_pattern,codon,BE_type):#Highest
    code={'UAU': 'Tyr', 'UGG': 'Trp', 'ACC': 'Thr', 'GAC': 'Asp', 'AGU': 'Ser', 'AGA': 'Arg', 'GGG': 'Gly', 'GCA': 'Ala', 'CGU': 'Arg', 'UCU': 'Ser', 'GGC': 'Gly', 'CCU': 'Pro', 'CAA': 'Gln', 'UGC': 'Cys', 'CAU': 'His', 'ACU': 'Thr', 'GUA': 'Val', 'UUG': 'Leu', 'AUG': 'Met', 'CUC': 'Leu', 'CAG': 'Gln', 'GUG': 'Val', 'UUC': 'Phe', 'UUU': 'Phe', 'UGU': 'Cys', 'CGG': 'Arg', 'GUC': 'Val', 'AGC': 'Ser', 'CCA': 'Pro', 'ACG': 'Thr', 'GAU': 'Asp', 'AAU': 'Asn', 'CGC': 'Arg', 'UAG': 'Stop', 'CGA': 'Arg', 'UAA': 'Stop', 'GCT': 'Ala', 'CCC': 'Pro', 'UGA': 'Stop', 'AUC': 'Ile', 'GAA': 'Glu', 'UCA': 'Ser', 'CUA': 'Leu', 'UCG': 'Ser', 'CUG': 'Leu', 'GUU': 'Val', 'CAC': 'His', 'GCG': 'Ala', 'AUA': 'Ile', 'AGG': 'Arg', 'GAG': 'Glu', 'AAA': 'Lys', 'UUA': 'Leu', 'AAC': 'Asn', 'AAG': 'Lys', 'GCU': 'Ala', 'CUU': 'Leu', 'UCC': 'Ser', 'UAC': 'Tyr', 'GGU': 'Gly', 'ACA': 'Thr', 'AUU': 'Ile', 'CCG': 'Pro', 'GGA': 'Gly'}
    start=seq_range.find(target_seq)
    wb4=pandas.ExcelWriter('Chart.xlsx')
    df=pandas.DataFrame(columns=('A','T','G','C'))
    for x in range(len(Substitution_pattern)):
        position_sub=0
        for i in Substitution_pattern[x]:
            position_sub+=i
        for y in range(4):
            wt_nt=int(seq_range[x].translate(seq_range[x].maketrans('ATGC','0123')))
            if wt_nt==y:
                Substitution_pattern[x][y]=round((all_c-position_sub)*100/all_c,2)
            else:
                Substitution_pattern[x][y]=round(Substitution_pattern[x][y]*100/all_c,2)
    for i in range(len(target_seq)):
        df.loc[i]=Substitution_pattern[i]
    df=df.T
    df.to_excel(wb4,'Sheet1',startrow=1)
    workbook=wb4.book
    ws4=wb4.sheets['Sheet1']
    cell_format=workbook.add_format({'align':'center','valign':'vcenter','bold':True})
    middle=workbook.add_format({'align':'center','valign':'vcenter'})
    blue=workbook.add_format({'align':'center','valign':'vcenter','bold':True,'font_color':'blue'})
    cor_site=workbook.add_format({'align':'center','valign':'vcenter','bg_color':'cyan'})
    BE_cell=workbook.add_format({'align':'center','valign':'vcenter','bold':True,'font_color':'red'})
    ws4.write_column(2,0,['A','T','G','C'],cell_format)
    for i in range(len(target_seq)):
        ii=i+start
        ws4.write(1,i+1,target_seq[i],cell_format)
    mut_target=''
    for i in range(len(target_seq)):
        ii=i+start
        seq_nt=Highest(Substitution_pattern[ii])
        mut_target+=str(seq_nt).translate(str(seq_nt).maketrans('0123','ATGC'))
        wt_nt=int(target_seq[i].translate(target_seq[i].maketrans('ATGC','0123')))
        ws4.write(seq_nt+2,i+1,Substitution_pattern[ii][seq_nt],cor_site)            
    for i in range(len(target_seq)-3,len(target_seq)):
        ws4.write(1,i+1,target_seq[i],blue)
    red=workbook.add_format({'color':'red'})
    green=workbook.add_format({'color':'green'})
    for i in range(8):
        if i*3+codon+2>len(target_seq):
            continue
        ws4.merge_range(0,i*3+codon+1,0,i*3+codon+3,'')
    #    if code[target_seq[i*3+codon:i*3+codon+3].replace('T','U')] != code[mut_target[i*3+codon:i*3+codon+3].replace('T','U')]:
    #        ws4.write_rich_string(0,i*3+codon+1,'',red,code[target_seq[i*3+codon:i*3+codon+3].replace('T','U')],u'→',green,code[mut_target[i*3+codon:i*3+codon+3].replace('T','U')],cell_format)
    #    else:
        ws4.write_rich_string(0,i*3+codon+1,code[target_seq[i*3+codon:i*3+codon+3].replace('T','U')],cell_format)
    if target_seq[4]=='C':
        ws4.write(1,5,'C',BE_cell)
        ws4.write_rich_string('F1','',red,code[target_seq[4:7].replace('T','U')],u'→',green,code['U'+mut_target[5:7].replace('T','U')],cell_format)
    wb4.save()
    wb4.close()

    return 0   

def Substitution_fold(end_range,seq_range,target_seq,data_c,con_c,data_list,control_list):
    wb5=pandas.ExcelWriter('Substitution_control_compare.xlsx')
    workbook=wb5.book
    df=pandas.DataFrame(columns=('control','mutation rate'))
    df1=pandas.DataFrame(columns=('data','mutation rate'))
    df2=pandas.DataFrame(columns=('fold',''))
    row=0
    start=seq_range.find(target_seq)
    for i in range((-1)*(start+17),len(seq_range)-(start+17)+1):
        if i==0:
            continue
        con_mut_count=0
        for x in range(4):
            con_mut_count+=control_list[row][x]
        con_rate=round(con_mut_count*100/con_c,2)
        data_rate=100-data_list[row][int(seq_range[row].translate(seq_range[row].maketrans('ATGC','0123')))]
        df.loc[i]=[seq_range[row],con_rate]
        df1.loc[i]=[seq_range[row],data_rate]
        con_mut_count=0
        for x in range(4):
            con_mut_count+=control_list[i][x]
        if con_rate != 0:
            df2.loc[i]=[round(data_rate/con_rate,2),'']
        elif con_rate==0:
            df2.loc[i]=[1.1,'']        
        row+=1
    df.T.to_excel(wb5,'Sheet1')
    df1.T.to_excel(wb5,'Sheet1',startrow=4)
    df2.T.to_excel(wb5,'Sheet1',startrow=8)
    workbook=wb5.book
    ws5=wb5.sheets['Sheet1']
    col=0
    ws5.write(10,0,'')
    chart=workbook.add_chart({'type':'scatter'})
    chart.add_series({'categories':['Sheet1',8,1,8,len(seq_range)],'values':['Sheet1',9,1,9,len(seq_range)]})
    chart.set_title({'name':'compare data and control'})
    chart.set_x_axis({'min':(-1)*(start+17),'max':len(seq_range)-(start+17)})
    ws5.insert_chart('A12',chart)
    wb5.save()
    wb5.close()

####### input
file_name='test.fastqjoin'
control_file='testc.fastqjoin'

wt_seq='aggggataccaccgatctctgtgatctgcgactgttttctctgtctgtgcaggtccacagtatggcattgcccgtgaagatgtggtcctgaatcgtattcttggggaaggcttttttggggaggtctatgaaggtgtctacacaaatcacgtgagttctaggatcttcccttacactcctcttccacatgtctgtagggtgagacagagctcgaa'.upper()
target_seq='GGTCCTGAATCGTATTCTTGGGG'.upper()

direc=os.getcwd()
codon = 1 # 0,1,2

filt_r=5
filt_n=1
end_range=70
indicator_range=15
BE_type=0

########

test_wt='cagtgaggatgccacgtagaaataggaggcctggccttgctcaaccacagttacacgggtgcccaatgggctggccacacagtcggggcagtcatcgggccggttatggtgggctgagaagaggcaggctggcgctgccagatgcacggctgtcccttggggctctaggtcatgcaggaagcagcggccctgcaggctggagccacaactgaccag'.upper()
reverse_wt='CTGGTCAGTTGTGGCTCCAGCCTGCAGGGCCGCTGCTTCCTGCATGACCTAGAGCCCCAAGGGACAGCCGTGCATCTGGCAGCGCCAGCCTGCCTCTTCTCAGCCCACCATAACCGGCCCGATGACTGCCCCGACTGTGTGGCCAGCCCATTGGGCACCCGTGTAACTGTGGTTGAGCAAGGCCAGGCCTCCTATTTCTACGTGGCATCCTCACTG'.upper()
complementary_target='GTCATCGGGCCGGTTATGGTGGG'.upper()
incorrect_target='CCCCCCCCCCCCCCCCCCCCCCC'.upper()
complementary_seq_range='TGCTCAACCACAGTTACACGGGTGCCCAATGGGCTGGCCACACAGTCGGGGCAGTCATCGGGCCGGTTATGGTGGGCTGAGAAGAGGCAGGCTGGCGCTGCCAGATGCACGGCTGTCCCTTGGGGCTCTAGGTCATGCAG'
reverse_seq_range='CTGCATGACCTAGAGCCCCAAGGGACAGCCGTGCATCTGGCAGCGCCAGCCTGCCTCTTCTCAGCCCACCATAACCGGCCCGATGACTGCCCCGACTGTGTGGCCAGCCCATTGGGCACCCGTGTAACTGTGGTTGAGCA'.upper()
complementary_s_seq='GTTATGGTGG'
reverse_s_seq='CCACCATAAC'

test_seq_range='ATGCGTA'
test_indicator_range=2
test_for_table=['AT','AT','TT','GT','CT','AA','AT','AG','AC']
test_back_table=['TA','AA','TA','GA','CA','TA','TT','TG','TC']

(test_for_table2,test_back_table2)=Indicator_table('ATGCCGTAACGTAGCTGA',6)
test_correct_seq='AAGCCGTAACGTAGCTAA'
test_for_seq='ATGAAGTAACGTAGCTGA'
test_back_seq='ATGCCGTAACGTAGAAGA'

emboss_wt='AAGCCGTAACGTAGCTAA'
sub_seq='AAGCCGTAAGGTAGCTAA'
del_seq='AAGCCGTAA-GTAGCTAA'
sub_seq='AAGCCGTAAgGTAGCTAA'
sub_sym='|||||||||.||||||||'
del_seq='AAGCCGTAA-GTAGCTAA'
del_sym='||||||||| ||||||||'

test_file='test.fastqjoin'
control_file='testc.fastqjoin'
test_seq_range2='ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC'
test_line_list={'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 4, 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTTTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 2, 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTATTAAAAAAACTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 4, 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTTTAAAGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 2, 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCAAATTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 2, 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCTTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 2, 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTAGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': 2}
test_emboss_dict={'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTTTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTtTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 2, 'Substituion'], 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 4, 'WT'], 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTATTAAAAAAACTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATT-------CTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||       ||||       |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCC-------TATTAAAAAAACTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 4, 'Others'], 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCAAATTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTC---TTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCAAATTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 2, 'Insertion'], 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCTTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCtTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 2, 'Substituion'], 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTTTAAAGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||| ||...||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGT-TTaaaGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 2, 'Deletion'], 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTAGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC': ['ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTATTCTTGGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', '|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||', 'ACTGTTTTCTCTGTCTGTGCAGGTCCACAGTATGGCATTGCCCGTGAAGATGTGGTCCTGAATCGTA-----GGGGAAGGCTTTTTTGGGGAGGTCTATGAAGGTGTCTACACAAATCATGTGAGTTCTAGGATCTTCCC', 2, 'Deletion']}
test_Count_list=[4, 4, 4, 2, 4, 4]
test_Substitution_position=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
test_Substitution_pattern=[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
test_CtoD=[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

def test_answer():

    #Seq_standard
    assert Seq_standard(test_wt,complementary_target,filt_r,end_range)==(complementary_seq_range,complementary_s_seq)
    assert Seq_standard(reverse_wt,complementary_target,filt_r,end_range)==(reverse_seq_range,reverse_s_seq)
    assert Seq_standard(test_wt,incorrect_target,filt_r,end_range)==(test_wt,'None')

    #Indicator_table
    assert Indicator_table(test_seq_range,test_indicator_range)==(test_for_table,test_back_table)

    #Seqcheck
    assert Seq_check(test_correct_seq,test_for_table2,test_back_table2)==test_correct_seq
    assert Seq_check(test_for_seq,test_for_table2,test_back_table2)==0
    assert Seq_check(test_back_seq,test_for_table2,test_back_table2)==0

    #substitution check
    assert Substitution(emboss_wt,sub_seq)==(emboss_wt,sub_seq,sub_sym)
    assert Substitution(emboss_wt,del_seq)==(emboss_wt,del_seq,del_sym)

    (t_for_table,t_back_table)=Indicator_table(test_seq_range2,15)
    #Read_seq
    assert Read_seq(test_file,t_for_table,t_back_table,filt_n)==(test_line_list, 22, 19, 1)

    emboss_needle(test_seq_range2, test_line_list, direc)
    #Read_EMBOSS
    assert Read_EMBOSS('emboss_after.txt',test_line_list,test_seq_range2)==(test_emboss_dict, test_Count_list, test_Substitution_position, test_Substitution_pattern, test_CtoD)

    (test_emboss_dict, test_Count_list, test_Substitution_position, test_Substitution_pattern, test_CtoD)==Read_EMBOSS('emboss_after.txt',test_line_list,test_seq_range2)
    #Total_Substitutions
    Total_Substitutions(test_seq_range2,target_seq,test_Substitution_position,test_Substitution_pattern)
    wb=pandas.read_excel('Total_Substitutions.xlsx').values.tolist()
    wbt=pandas.read_excel('test_Total_Substitutions.xlsx').values.tolist()
    for i in range(len(wb)):
        wb[i].remove(wb[i][3])
        wbt[i].remove(wbt[i][3])
    assert wb==wbt

    #C_to_D_substitution
    C_to_D_Substitutions(test_seq_range2,target_seq,test_CtoD)
    wb=pandas.read_excel('C_to_D_Substitutions.xlsx').values.tolist()
    wbt=pandas.read_excel('test_C_to_D_Substitutions.xlsx').values.tolist()
    for i in range(len(wb)):
        wb[i].remove(wb[i][3])
        wbt[i].remove(wbt[i][3])
    assert wb==wbt

    #Transition_figure
    Transition_figure(test_seq_range2,test_Count_list[1],test_Substitution_pattern)
    wb=pandas.read_excel('Substitution_transition.xlsx').values.tolist()
    wbt=pandas.read_excel('test_Substitution_transition.xlsx').values.tolist()
    assert wb==wbt

    #Highest
    assert Highest([1,5,2,40])==3

    #Chart
    Chart(test_seq_range2,target_seq,18,test_Substitution_pattern,codon,BE_type)
    wb=pandas.read_excel('Chart.xlsx').values.tolist()
    wbt=pandas.read_excel('test_Chart.xlsx').values.tolist()
    assert wb==wbt

    (control_line_list,control_all_count,control_primer_count,control_filted_count)=Read_seq(control_file,t_for_table,t_back_table,filt_n)
    emboss_needle(test_seq_range2, control_line_list,direc)
    (control_emboss_dict,control_Count_list,control_Substitution_position,control_Substitution_pattern,CtoD)=Read_EMBOSS('emboss_after.txt', control_line_list, test_seq_range2)
    #Substitution_fold
    Substitution_fold(end_range,test_seq_range2,target_seq,18,21,test_Substitution_pattern,control_Substitution_pattern)
    wb=pandas.read_excel('Substitution_control_compare.xlsx').values.tolist()
    wbt=pandas.read_excel('test_Substitution_control_compare.xlsx').values.tolist()
    assert wb==wbt
########
