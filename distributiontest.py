# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 19:06:54 2020

@author: s4142554
"""

import os


from scipy.stats import shapiro, normaltest, anderson
import scipy.stats as stats
###############################################################################
script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
###############################################################################

def SkewTest(beta_list):

    skew = stats.skew(beta_list) if len(beta_list) >= 20 else "*"
    
    if skew == "*":
        skewness = "*"
    else:
        if skew >= -0.5 and skew <= 0.5:
            skewness = "symmetric" #changed from normal
        elif skew >= -1 and skew < -0.5:
            skewness = "moderate positive skew"
        elif skew > 0.5 and skew <= 1:
            skewness = "moderate negative skew"
        elif skew < -1:
            skewness = "high positive skew"
        elif skew > 1:
            skewness = "high negative skew"
    
    return skew, skewness

def KurtosisTest(beta_list):
    
    kValue = stats.kurtosis(beta_list) if len(beta_list) >= 20 else "*"#excess kurtosis

    if kValue == "*":
        shape = "*"
    else:
        if kValue > 2:
            shape = "Leptokuric"
        elif kValue <=2 and kValue >= -2:
            shape = "Mesokuric" #normal dist
        elif kValue < -2:
            shape = "Platykuric"
    return kValue, shape

def newNormalTest(dataList):
    
    if len(dataList) >= 10: #sample saize must be more than or equal to 10 to proceed
        #1 shapiro-wilk test
        shapiro_stat, shapiro_p = shapiro(dataList)
        
        #2 D'Agostino's K^2 Test
        dagostino_stat, dagostino_p = normaltest(dataList)
        
        #3 Skewness
        skew, skewness = SkewTest(dataList)
        
        #4 Kurtosis
        kVal, kurtosis = KurtosisTest(dataList)
        
        #5 Anderson-Darling test
        #if stat < cv - data is normal
        #if stat > cv - data follows non normal dis
        and_res = anderson(dataList)
        anderson_stat = and_res.statistic
        anderson_cv = and_res.critical_values #list of critical values
        sig_level = and_res.significance_level

        andList = [anderson_stat]

        #go through list of critical values and find significance at each level
        for i in range(len(anderson_cv)):
            sl, cv = sig_level[i], anderson_cv[i]
            
            if anderson_stat <= anderson_cv[i]:
    #            norm =  sl, cv, "normal"
    #            andList.extend(norm) #sl, cv
                andList.extend([cv, "normal"])
                
            else:
    #            norm = sl,cv, "non-normal"
    #            andList.extend(norm) #sl, cv
                andList.extend([cv, "non-normal"])                 

    else:
        shapiro_stat = "NA"
        shapiro_p = "NA"
        dagostino_stat = "NA"
        dagostino_p = "NA"
        
        andList = ["NA" for n in range(11)]
        
    result = [shapiro_stat, shapiro_p, dagostino_stat, dagostino_p, skew, skewness, kVal, kurtosis]
    result.extend(andList)
    
    #result = [shapiro_stat, shapiro_p, dagostino_stat, dagostino_p, skew, skewness, kVal, kurtosis, 
    #anderson_stat, anderson_cv15, anderson_dis15, anderson_cv10, anderson_dis10,
    #anderson_cv5, anderson_dis5, anderson_cv2.5, anderson_dis2.5, anderson_cv1, anderson_dis1]

    
    return result # results report the stats of shapiro-wilk, dagostino, skewness, kurtois, and anderson-darling



def FDR(pVal_List, method, vals=None):
    
    #uses the benjamini hochberg method
    #revised FDR to report adjusted FDR, new version takes into account NA values
    if method == "BH":
#        sorted_pVal_List = [i for i in sorted((pVal_List), key=lambda x:x[1])]
#        sorted_pVal_List = [i for i in sorted((pVal_List), key=lambda x: x[1])]
        sorted_pVal_List =  [i for i in sorted(enumerate(pVal_List), key=lambda x:x[1])]

        ori_p_val_index = [i[0] for i in sorted_pVal_List] #orginial index of sorted list
        sorted_p_values = [i[1] for i in sorted_pVal_List] #get sorted pVal list
        
        # 2 rank each p values
        p_value_ranks = [n+1 for n in range(len(sorted_p_values))] #rank each one according to place in list
        
        #3 largest adjust p value == largest p value
        p_adj = sorted_p_values[-1]

        adj_p_values = [[ori_p_val_index[-1], p_adj, p_adj]] #[index, p val, adjusted p val]

        
        #4 find adjusted p values here, start from the second p val as the first is already in the adj_p_values list
        for ind, p, r in zip(ori_p_val_index[::-1][1:], sorted_p_values[::-1][1:], p_value_ranks[::-1][1:]):
            
            

            p_adj = p*(len(sorted_p_values)/r) # current adjusted p val =  p val * (total # of p vals/ p value rank)

            #new adjusted p val = min(p_adj, find min p_adj between p_adj and previous p_adj
            #previous p is the last item in the adj_p_values
#                    new_p_adj = round(min(p_adj, adj_p_values[-1]), 2)
            new_p_adj = round(min(p_adj, adj_p_values[-1][2]), 2)

            #add to adj_p_values
            adj_p_values.append([ind, p, new_p_adj]) #[ori_i, p_Val, adjust_pVal]


    
    for line in adj_p_values:
        ind, pval, p_adj_val = line[0], line[1], line[2]
        if pval <= p_adj_val: 
            pval, p_adj_val = pval, p_adj_val
            break

    #resort pVals by the original index
    sorted_adj_pVals = [i for i in sorted((adj_p_values), key=lambda x:x[0])]
#    result = [pval, adj_p_values, sorted_adj_pVals]

    if vals == None:
        result = [pval, p_adj_val]
        
    else:
        result = [pval, adj_p_values, sorted_adj_pVals]

    return result


def HypothesisTest(p_val, alpha):
    if p_val == "NA":
        result = "NA"
    elif alpha == -1:
        result = "-1"
    else:
        p_val = float(p_val)
        alpha = float(alpha)

        if p_val > alpha: #accept null hypothesis
            result = "insignificant"#"normal" 
        elif p_val <= alpha: #reject null hypothesis
            result = "significant" #"non-normal"

    return result



def DistScore(skewness, kurtosis, shapiro, dagostino, andersonDist):
    
    DistDict = {
            0: "non-Gaussian",
            1: "non-Gaussian-like",
            2: "non-Gaussian-like",
            3: "Gaussian-like",
            4: "Gaussian-like",
            5: "Gaussian"
            }
    
    dScore = 0
    #skewness, kurtosis, shapiro, dagostino, andersonDist
    dScore+=1 if skewness == "symmetric" else 0
    dScore+=1 if kurtosis == "Mesokuric" else 0
    dScore+=1 if shapiro == "normal" else 0
    dScore+=1 if dagostino == "normal" else 0
    dScore+=1 if andersonDist == "normal" else 0
    
    
    distribution = DistDict[dScore]
    
    return dScore, distribution



#uses the nromality test to predict the distribution of the matrix
def MatrixDistributionTest(InputMatrixFile, ResultFile): 
    
    tallyDict = {
        "dScore": ["distribution", "tally"],
            0: ["non-Gaussian", 0],
            1: ["non-Gaussian-like", 0],
            2: ["non-Gaussian-like", 0],
            3: ["Gaussian-like", 0],
            4: ["Gaussian-like", 0],
            5: ["Gaussian", 0]
        }

    title = [
            "ProbeID", "n", "shapiro_stat", "shapiro_p", "dagostino_stat", "dagostino_p", "skew", "skewness", "kVal", "kurtosis",
            "anderson_stat", "anderson_cv15", "anderson_dis15", "anderson_cv10", "anderson_dis10",
            "anderson_cv5", "anderson_dis5", "anderson_cv2.5", "anderson_dis2.5", "anderson_cv1", "anderson_dis1", 
            ]
    
    ncol = len(title) - 2

    #result = [shapiro_stat, shapiro_p, dagostino_stat, dagostino_p, skew, skewness, kVal, kurtosis, 
    #anderson_stat, anderson_cv15, anderson_dis15, anderson_cv10, anderson_dis10,
    #anderson_cv5, anderson_dis5, anderson_cv2.5, anderson_dis2.5, anderson_cv1, anderson_dis1]    

    name = InputMatrixFile.split('\\')[-1]
    
   
    #create temperory file in cwd
    temp_file =  "tempFile.txt"
    
    
    print ("Processing...", name)
    
    ncol = len(title) - 2

    Shap_pValues = [] #empty list for adding p vals
    dagostino_pValues = [] #empty list for adding p vals
        
    with open(InputMatrixFile, 'r') as mfile, open(temp_file, 'w') as tfile:
        for i, line in enumerate(mfile):
            line = line.split()
            probe = line[0]
            bVal = line[1:]
            
            # start processing from the second line of the file. (ie. samples from the first probe)
            if i >=1:
                #proceed with checking
#                print (i, probe)
                
                if len(set(bVal)) == 1:
                    rsline = [probe, 0]
                    rsline.extend(['NA' for n in range(ncol)])
                elif len([float(n) for n in bVal if n != "NA"]) < 9: #if contents less than 10
                    rsline = [probe, len([float(n) for n in bVal if n != "NA"])]
                    rsline.extend(['NA' for n in range(ncol)])       
                else:
                    bVal = [float(n) for n in bVal if n != "NA"]
                    norm = newNormalTest(bVal) 
                    
                    
                    #add all pValues to list for p value cut off check
                    Shap_pValues.append(norm[1])
                    dagostino_pValues.append(norm[3])

                    
                    rsline = [probe, len(bVal)]
                    rsline.extend(norm)
                
                rsline = '\t'.join([str(n) for n in rsline])+'\n'
                tfile.write(rsline) #write results to temp file 


    title = [
            "ProbeID", "n", "shapiro_stat", "shapiro_p", "dagostino_stat", "dagostino_p", "skew", "skewness", "kVal", "kurtosis",
            "anderson_stat", "anderson_cv15", "anderson_dis15", "anderson_cv10", "anderson_dis10",
            "anderson_cv5", "anderson_dis5", "anderson_cv2.5", "anderson_dis2.5", "anderson_cv1", "anderson_dis1", 
            "shapiro_dist", "dagostino_dist", "andersonDarl_dist", "dScore", "Distribution"
             ]
    
    title = "\t".join(title)+"\n"
                
    print ("Now determining p cutoff...")
    #find false discovvery rate 
    Shap_p_cutoff, Shap_adj_p_cutoff = FDR(Shap_pValues, "BH") #[pval, adj_p_values, sorted_adj_pVals]
    dagostino_p_cutoff, dagostino_adj_p_cutoff= FDR(dagostino_pValues, "BH") #[pval, adj_p_values]
    print ("Alpha values:", "Shapiro-Wilk:", Shap_adj_p_cutoff, "D'Agostino:", dagostino_p_cutoff)##
    print ("Calculating distributons...")

    with open(temp_file, 'r') as tfile, open(ResultFile, 'w') as rfile:
        rfile.write(title)
        for i, line in enumerate(tfile):
        
            line = line.split()
            
            probe = line[0]

            
#            print ("Calculating distributons", i, probe)
            
            shap_p = float(line[3]) if line[3] != "NA" else "NA"
            dagostino_p = float(line[5]) if line[5] != "NA" else "NA"
            
            skewness = line[7] if line[7] != "NA" else "NA"
            kurtosis = line[9] if line[9] != "NA" else "NA"

            anderson = [line[12], line[14], line[16], line[18], line[20]]
            
            if shap_p != "NA" and dagostino_p != "NA":
                shapiro = "normal" if HypothesisTest(shap_p, Shap_p_cutoff) == "insignificant" else "non-normal" 
                dagostino = "normal" if HypothesisTest(dagostino_p, dagostino_p_cutoff) == "insignificant" else "non-normal"
            
                andersonDist = "normal" if anderson.count("normal") in [4,5] else "non-normal"

                dScore, distribution = DistScore(skewness, kurtosis, shapiro, dagostino, andersonDist)
                
                #add dScore to tally
                tallyDict[dScore][1]+=1
                
                rline = line
                rline.extend([andersonDist, dScore, distribution])
                
            else:
                rline = line
                rline.extend(["NA" for n in range(5)])
                
            rline = "\t".join([str(n) for n in rline])+"\n"
            rfile.write(rline)
            
    #remove temp file
    os.remove(temp_file)
    
    print ("Finished processing...\n")
    #print summary
    for key, value in tallyDict.items():
        print (key, value[0], value[1])
    return ""