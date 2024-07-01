from array import array

jetEtaBins = [0., 1.3, 2.5, 3., 3.5, 4., 5.]
egEtaBins = [0., 1.479, 2.5]
muEtaBins = [0., 0.83, 1.24, 2.4]
muEMTFBins = [1.24, 1.6, 2.1, 2.4]


ht_bins = array('f', [ i*10 for i in range(50) ] + [ 500+ i*20 for i in range(25) ] + [1000 + i*50 for i in range(10)] +[1500,1600,1700,1800,2000,2500,3000,4000,6000,8000,14000])
#leptonpt_bins = array('f',[ i for i in range(50) ] + [ 50+2*i for i in range(10) ] + [ 70+3*i for i in range(10) ] + [100+10*i for i in range(10) ] + [200, 250, 300, 400, 500])
jetmetpt_bins = array('f',[ i*5 for i in range(50) ] +  [250+10*i for i in range(25) ]  + [500+20*i for i in range(10) ] + [700, 800, 900, 1000, 1200, 1500, 2000, 3000 ])

# Finer binning:
leptonpt_bins = array('f',[i * .2 for i in range(15 * 5)] + [ 15 + i * .5 for i in range(5 * 2)] + [20 + i for i in range(30)] + [ 50+2*i for i in range(10) ] + [ 70+3*i for i in range(10) ] + [100+10*i for i in range(10) ] + [200, 250, 300, 400, 500])

coarse_leptonpt_bins = array('f',[ i for i in range(50) ] + [ 50+2*i for i in range(10) ] + [ 70+3*i for i in range(10) ] + [100+10*i for i in range(10) ] + [200, 250, 300, 400, 500, 700, 1000])
coarse2_leptonpt_bins = array('f',[ i for i in range(20) ] + [20, 25, 30, 35, 45, 60, 75, 100, 140, 200, 250, 300, 400, 500])

