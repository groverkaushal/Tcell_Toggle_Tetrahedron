import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr
import statistics
import os
from barplot_annotate import barplot_annotate_brackets




input_gseid = "GSE135390"


##### location of output directory
if not os.path.exists("./results/" + input_gseid):
    os.mkdir("./results/" + input_gseid)


#### ssgsea dataframe is sorted
data = pd.read_csv("./ssGSEAresults/gseapy.gene_set.ssgsea.report.csv")
data = data.sort_values(by=['Name','Term'])
data = data.reset_index()
data = data.iloc[:,1:5]



###### Print the DataFrame without automatically shortening it in console. 
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.expand_frame_repr', False)  # Do not wrap DataFrame in the console
pd.set_option('display.max_rows', None)  # Display all rows
pd.set_option('display.width', None)  # Display without line-wrapping




coordinates = []

for signature in ['TH1','TH2','TH17','TREG','TH']:


    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    temp5 = []
    temp6 = []
    temp7 = []
    # temp8 = []
    # temp9 = []
    # temp10 = []
    # temp11 = []
    # temp12 = []
#     # # # # temp13 = []
#     # # # # temp14 = []
#     # # # # temp15 = []
#     # # # # temp16 = []


    for i in range(len(data)):
        # temp = data.iloc[i,0][:3]
        # time = data.iloc[i,0][:2]
        # gfp_all = data.iloc[i,0][::-1][2:5][::-1]


        # if data.iloc[i,1]== signature and data.iloc[i,0][:10] == "Day0 Naive":
        #     temp1.append(data.iloc[i,3])
        
        # if data.iloc[i,1]== signature and data.iloc[i,0][::-1][:5][::-1] == "Naive":
        #     temp1.append(data.iloc[i,3])
        
        if data.iloc[i,1]== signature and data.iloc[i,0][::-1][:3][::-1] == "Th1":
            temp2.append(data.iloc[i,3])

        if data.iloc[i,1]== signature and  data.iloc[i,0][::-1][:6][::-1] == "Th1/17":
            temp7.append(data.iloc[i,3])

        if data.iloc[i,1]== signature and  data.iloc[i,0][::-1][:4][::-1] == "Th17":
            temp4.append(data.iloc[i,3])

        if data.iloc[i,1]== signature and data.iloc[i,0][::-1][:3][::-1] == "Th2":
            temp3.append(data.iloc[i,3])
        
        # if data.iloc[i,1]== signature and data.iloc[i,0][::-1][:4][::-1] == "Th22":
        #     temp6.append(data.iloc[i,3])

        if data.iloc[i,1]== signature and data.iloc[i,0][::-1][:5][::-1] in ['Treg1','g1/17','reg17','Treg2','reg22']:
            temp5.append(data.iloc[i,3])

        # if data.iloc[i,1]== signature and data.iloc[i,0][::-1][2:][::-1] == 'ExVivoFoxp3-RORÎ³t+':
        #     temp3.append(data.iloc[i,3])

        # if data.iloc[i,1]== signature and data.iloc[i,0][::-1][2:][::-1] == 'ExVivoFoxp3+RORÎ³t+':
        #     temp4.append(data.iloc[i,3])

        # if data.iloc[i,1]== signature and data.iloc[i,0][:6] == 'Th1_WT':
        #     temp5.append(data.iloc[i,3])
        
        # if data.iloc[i,1]== signature and data.iloc[i,0][:6] == 'Th1_KO':
  

    stddeviation = []
    # stddeviation.append(statistics.stdev(temp1))
    stddeviation.append(statistics.stdev(temp2))
    stddeviation.append(statistics.stdev(temp3))
    stddeviation.append(statistics.stdev(temp4))
    stddeviation.append(statistics.stdev(temp5))
    # stddeviation.append(statistics.stdev(temp6))
    stddeviation.append(statistics.stdev(temp7))
    # stddeviation.append(statistics.stdev(temp8))
    # stddeviation.append(statistics.stdev(temp9))
    # stddeviation.append(statistics.stdev(temp10))
    # stddeviation.append(statistics.stdev(temp11))
    # stddeviation.append(statistics.stdev(temp12))
#     # # # # stddeviation.append(statistics.stdev(temp13))
#     # # # # stddeviation.append(statistics.stdev(temp14))
#     # # # # stddeviation.append(statistics.stdev(temp15))
#     # # # # stddeviation.append(statistics.stdev(temp16))




    y_axis = []
    # y_axis.append(statistics.mean(temp1))
    y_axis.append(statistics.mean(temp2))
    y_axis.append(statistics.mean(temp3))
    y_axis.append(statistics.mean(temp4))
    y_axis.append(statistics.mean(temp5))
    # y_axis.append(statistics.mean(temp6))
    y_axis.append(statistics.mean(temp7))
    # y_axis.append(statistics.mean(temp8))
    # y_axis.append(statistics.mean(temp9))
    # y_axis.append(statistics.mean(temp10))
    # y_axis.append(statistics.mean(temp11))
    # y_axis.append(statistics.mean(temp12))
#     # # # # y_axis.append(statistics.mean(temp13))
#     # # # # y_axis.append(statistics.mean(temp14))
#     # # # # y_axis.append(statistics.mean(temp15))
#     # # # # y_axis.append(statistics.mean(temp16))

    xaxis = []
    for i in range(len(y_axis)):
        xaxis.append(i)



############# Plotting bar graph with GSM samples on xaxis and their signatures on y axis

    light_yellow = (255/255,255/255,153/255)
    light_green = (102/255, 255/255, 102/255)

    # colors = ['green', 'blue', 'pink', 'yellow', 'orange', 'red', 'purple']
    # labels = ['Naive', 'Th1', 'Th1/17', 'Th17', 'Th2', 'Th22', 'Treg']
    colors = colors = ['#ff856b','#65c8b4','#9ac6fa','#ffc700','#dda0fa']
    labels = ['Th1', 'Th17', 'Th2', 'Treg', 'Th1/17']


    from scipy.stats import ttest_ind
    import numpy as np

    pvalue_temp = []
    # pvalue_temp.append(temp1)
    pvalue_temp.append(temp2)
    pvalue_temp.append(temp3)
    pvalue_temp.append(temp4)
    pvalue_temp.append(temp5)
    # pvalue_temp.append(temp6)
    pvalue_temp.append(temp7)
    # pvalue_temp.append(temp8)
    # pvalue_temp.append(temp9)
    # pvalue_temp.append(temp10)
    # pvalue_temp.append(temp11)
    # pvalue_temp.append(temp12)



###### initialize an empty len(pvalue_temp)*len(pvalue_temp) dataframe 
    pvaluematrix = np.zeros((len(pvalue_temp),len(pvalue_temp)))


###### fill it with the pvalues taking 2 columns each time
    for i in range(len(pvalue_temp)):
        for j in range(len(pvalue_temp)):
            corr_coef, p_value = ttest_ind(pvalue_temp[i], pvalue_temp[j], equal_var=False)
            pvaluematrix[i,j] = p_value


    pvaluematrix = pd.DataFrame(pvaluematrix)
    print(pvaluematrix)


###### create a list which has tuples which contain the column number if the pvalue is less than 0.05
    indices_star = pvaluematrix.where((pvaluematrix > 0.01) & (pvaluematrix < 0.05)).stack().index
    indices_star = list(set(tuple(sorted(pair)) for pair in indices_star))

    indices_doublestar = pvaluematrix.where((pvaluematrix > 0.001) & (pvaluematrix < 0.01)).stack().index
    indices_doublestar = list(set(tuple(sorted(pair)) for pair in indices_doublestar))

    indices_triplestar = pvaluematrix.where(pvaluematrix < 0.001).stack().index
    indices_triplestar = list(set(tuple(sorted(pair)) for pair in indices_triplestar))

    # print(indices_star)



    plt.bar(xaxis, y_axis, yerr=stddeviation, capsize=4, color = colors)



##### barplot_annotate_brackets is imported into this code to plot lines between the columns and displaying asterisks for significant columns.
    tttt = 0.03
    for singlestar in indices_star:
        barplot_annotate_brackets(singlestar[0], singlestar[1], '*', xaxis, y_axis, dh=tttt, barh=.01)
        tttt = tttt + 0.078


    for doublestar in indices_doublestar:
        barplot_annotate_brackets(doublestar[0], doublestar[1], '**', xaxis, y_axis, dh=tttt, barh=.01)
        tttt = tttt+ 0.078


    for triplestar in indices_triplestar:
        barplot_annotate_brackets(triplestar[0], triplestar[1], '***', xaxis, y_axis, dh=tttt, barh=.01)
        tttt = tttt+0.078

    plt.xticks(xaxis, labels)
    plt.ylabel(signature)
    plt.title(input_gseid + " " + signature + " signature")
####### to display the mean value of the column on the each bar plots 



    patches = []
    for color, label in zip(colors, labels):
        patch = mpatches.Patch(color=color, label=label)
        patches.append(patch)

    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.yticks(np.linspace(-1.0, 0.7, 18))

    aspect_ratio = 1.5  # Width to height ratio
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(10, 13 / aspect_ratio)  # Set the figure size
    fig.tight_layout(rect=(0, 0, 0.98, 1))
    plt.savefig("./results/" + input_gseid + "/1_1." + signature + ".png")
    # plt.show()
    plt.clf() 






########### code to plot heatmap using the pvaluematrix 
    
    from matplotlib.colors import LinearSegmentedColormap
    colors_1 = [(0, 'red'), (0.05, 'grey'), (1, 'grey')]
    cmap = LinearSegmentedColormap.from_list('CustomMap', colors_1)

    # Plot the p-value heatmap with annotations using the 'coolwarm' colormap
    plt.imshow(pvaluematrix, cmap=cmap, interpolation='nearest')

    # Add annotations to each cell
    for i in range(pvaluematrix.shape[0]):
        for j in range(pvaluematrix.shape[1]):
            plt.text(j, i, f'{pvaluematrix.iloc[i, j]:.2f}', ha='center', va='center', color='w')

    plt.colorbar()
    plt.title(input_gseid + " " + signature + ' p-value Heatmap')

    # Set the custom labels for x-axis and y-axis
    plt.xticks(range(pvaluematrix.shape[1]), labels)
    plt.yticks(range(pvaluematrix.shape[0]), labels)

    aspect_ratio = 1  # Width to height ratio
    fig = plt.gcf()  # Get the current figure
    fig.set_size_inches(10, 10 / aspect_ratio)  # Set the figure size
    plt.savefig("./results/" + input_gseid + "/1_1." + signature + "_heatmap.png")
    # plt.show()
    plt.clf() 



    coordinates.append(pvalue_temp)






############ code to plot the scatter plot of all the GSM samples with ssgsea values of the one signature on xaxis and the ssgsea values of another signature on yaxis
    
from itertools import product
colors_patch = colors.copy()
labels_patch = labels.copy()
nested_lengths = [len(sublist) for sublist in coordinates[3]]
colors = [color for color, count in zip(colors, nested_lengths) for _ in range(count)]
labels = [color for color, count in zip(labels, nested_lengths) for _ in range(count)]


sigg = [0,1,2,3,4]
sigg_name = ['TH1','TH2','TH17','TREG','TH']
pairs = list(product(sigg, repeat=2))

for i in pairs:

    tt_temp_1 = i[0]
    tt_temp_2 = i[1]

    x_term = [item for sublist in coordinates[tt_temp_1] for item in sublist]
    y_term = [item for sublist in coordinates[tt_temp_2] for item in sublist]

    # Calculate linear regression line
    slope, intercept = np.polyfit(x_term, y_term, 1)

    # Generate more points for the regression line
    x_regression = np.linspace(min(x_term), max(x_term), 1000)
    y_regression = slope * x_regression + intercept


    # from scipy.stats import ttest_ind
    from scipy.stats import pearsonr
    corr_coef, p_value = pearsonr(x_term, y_term)


    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = fig.add_subplot(111)

    ax.scatter(x_term, y_term, color = colors,s=90)

    ax.set_xlabel(sigg_name[tt_temp_1])
    ax.set_ylabel(sigg_name[tt_temp_2])
    ax.set_title(input_gseid)
    ax.plot(x_regression, y_regression, color='red', label='Linear Regression Line')


    # Add the correlation coefficient value to the plot
    correlation_text = f"R : {corr_coef:.2f}"

    ax.text(0.89, 0.97, correlation_text, transform=ax.transAxes, fontsize=10)

    if p_value < 0.0001:
        p_value_text = "< 10^-4"
        correlation_text_2 = f"pvalue {p_value_text}"
        ax.text(0.80, 0.94, correlation_text_2, transform=ax.transAxes, fontsize=10)

    elif 0.0001<p_value<0.001:
        p_value_text = "< 10^-3"
        correlation_text_2 = f"pvalue {p_value_text}"
        ax.text(0.80, 0.94, correlation_text_2, transform=ax.transAxes, fontsize=10)

    elif 0.001<p_value<0.05:
        p_value_text = "< 0.05"
        correlation_text_2 = f"pvalue {p_value_text}"
        ax.text(0.82, 0.94, correlation_text_2, transform=ax.transAxes, fontsize=10)

    else:
        p_value_text = str(f"= {p_value:.2f}")
        correlation_text_2 = f"pvalue {p_value_text}"
        ax.text(0.82, 0.94, correlation_text_2, transform=ax.transAxes, fontsize=10)



    patches = []
    for color, label in zip(colors_patch, labels_patch):
        patch = mpatches.Patch(color=color, label=label)
        patches.append(patch)

    plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc='upper left',fontsize = 8)
    plt.yticks(np.linspace(-0.9, 0.1, 11))
    aspect_ratio = 1.2  # Width to height ratio
    fig.set_size_inches(8, 8 / aspect_ratio)  # Set the figure size
    fig.tight_layout(rect=(0, 0, 0.98, 1))

    plt.savefig("./results/" + input_gseid + "/_"+ sigg_name[tt_temp_1] + "vs" + sigg_name[tt_temp_2] +".png")
    # plt.show()
    plt.clf()




###### just printing all the data ultimately used for plotting all the graphs
coordinates = pd.DataFrame(coordinates,columns=labels_patch, index=['TH1','TH2','TH17','TREG','TH'])
coordinates.to_csv("./results/" + input_gseid + "/coordinates.txt",sep="\t")