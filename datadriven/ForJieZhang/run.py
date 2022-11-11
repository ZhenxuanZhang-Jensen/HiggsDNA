def decimal_range(start, stop, increment):
    while start < stop: # and not math.isclose(start, stop): Py>3.5
        yield start
        start += increment
import os
for i in decimal_range(0.1, 0.9, 0.1):
    # print("python smooth_DNN.py -d ./ --min "+str(i)+" --max 1. -n 30 -r 115 -R 135 --outDir smooth_DNN_%s/"%(str(i)))
    os.system("python smooth_DNN_WWgg.py -d ./ --min "+str(i)+" --max 1. -n 30 -r 115 -R 135 --outDir smooth_DNN_WWgg_%s/"%(str(i)))
    # for j in range(2,9):
        # print(i)
        # print(type(i))
    os.system("python opt.py -d ./ -n 30 -c 2 --massMin 115 --massMax 135 -s 1 -w 1 -m %s -o opt_WWgg"%(str(i)))
    os.system("mv opt_WWgg/categorize_nBins_30_nCat_2_massMin115.0_massMax135.0_v2.txt opt_WWgg/categorize_nBins_30_nCat_2_massMin115.0_massMax135.0_v2_m%s.txt"%(str(i)))
        # os.system("python opt.py -d ./ -n 30 -c %i --massMin 115 --massMax 135 -s 1 -w 1 -m %s -o opt_%s"%(j,str(i),str(i)))
        # print("python opt.py -d ./ -n 30 -c %i --massMin 115 --massMax 135 -s 1 -w 1 -m %f -o opt_%s"%(j,str(i),str(i)))
    # sys.exit()

# python smooth_DNN.py -d ./ --min 0.2 --max 1. -n 30 -r 115 -R 135 --outDir smooth_DNN_m02/
# python opt.py -d ./ -n 30 -c 3 --massMin 115 --massMax 135 -s 1 -w 1 -m 0.2 -o opt_m02