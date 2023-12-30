import matplotlib.pyplot as plt
import numpy as np

def barplot_annotate_brackets(num1, num2, data, center, height, yerr=None, dh=.05, barh=.1, fs=20, maxasterix=None):
    if type(data) is str:
        text = data
    # else:
    #     text = ''
    #     p = .05

    #     while data < p:
    #         text += '*'
    #         p /= 10.

    #         if maxasterix and len(text) == maxasterix:
    #             break

    #     if len(text) == 0:
    #         text = 'n. s.'

    lx, ly = center[num1], 0.05
    rx, ry = center[num2], 0.05

    if yerr:
        ly += yerr[num1]
        ry += yerr[num2]

    # ax_y0, ax_y1 = plt.gca().get_ylim()
    # dh *= (ax_y1 - ax_y0)
    # barh *= (ax_y1 - ax_y0)

    if max(height) > 0 :
        y = max(height) + dh
    else:
        y =  dh
    

    barx = [lx, lx, rx, rx]
    bary = [y, y+barh, y+barh, y]
    mid = ((lx+rx)/2, y+barh-0.02)

    plt.plot(barx, bary, c='black', linewidth=0.8)

    kwargs = dict(ha='center', va='bottom')
    if fs is not None:
        kwargs['fontsize'] = fs

    plt.text(*mid, text, **kwargs,color='red')






# heights = [1.8, 2, 3,4,6]
# bars = np.arange(len(heights))

# plt.figure()
# plt.bar(bars, heights, align='center')
# # plt.ylim(0, 5)
# barplot_annotate_brackets(0, 1, .1, bars, heights)
# barplot_annotate_brackets(1, 2, .001, bars, heights)
# barplot_annotate_brackets(0, 2, 'p < 0.0075', bars, heights, dh=.2)
# plt.show()