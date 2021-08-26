#  © Shahram Talei @ 2019 The University of Alabama - All rights reserved.
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.


import matplotlib.pyplot as plt
import numpy as np
from matplotlib.legend_handler import HandlerLine2D
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D
import argparse

plt.rcParams["font.size"] =12

h_title='m12b,z=0'
#MWSat=np.array([132,108,66,42,47,209,218,160,105,36,28,116,30,251,183,44,117,24,76,22,380,147,120,53,132,79,78,151,28,214,100,145,50,254,233,154,178,215,83,114,46,182,30,92,26,69,84,23,35,86,62,30,58,25,48,55,97,32,76,91,38]) # Drlica-Wagner_2020
MWSat=np.array([24,28,29,23,19,26,29,35,32,33,41,37,43,43,40,46,50,46,45,50,48,52,54,61,64,63,78,79,82,81,86,87,87,87,92,98,102,94,108,105,116,116,128,133,126,141,148,155,161,170,182,188,203]) #Putman_2021
MWSat_sigma=np.array([2.9,3.7,1.2,11.4,2.7,5.6,3.4,3.3,5.6,3.4,3.4,4.6,10.5,4,10.7,4.3,1.9,8.6,27.6,4.6,9.5,9.1,9.2,4.9,7.9,7,6.6,5.4,2.7,2.9,5.7,5.1,11.7,3.6,3.3,4.6,2.3,5.4,5.4]) #Putman_2021
M31Sat=np.array([184,175,181,130,139,144,66,187,108,128,55,73,102,39,126,109,140,161,27,90,115,46,77,135,134,182]) #Putman_2021
#MWSat_m=np.array([19,0.91,6.3,0.94,11,56,2.6,138000,12,4.6,1.3,1.1,190,14,0.26,0.23,25,6500,11,3.9,9.5,0.27])#*10^6# M_sun
MWSat_m=np.array([74338.6888757733,152855.119951337,24787.3167488654,160947062.324042,179749.916896755,437685.953403389,86048.4629480132,258382.654826433,2027944.9174357,494778.661951076,215121.157370033,679250.637124876,2000361.58391393,245640.075889657,1092460.97614843,1438622739.9121,156190.800527338,4163394.58626644,391602758.62069,1880244.51725872,17004643.6057989,8899386.42443484,10986140.739585,446804.410765959,13239888.4358726,6724862.30489582,6302510.24713887,2170565.03422497,3615347.38513096,109551.752027455,43819077.1672987,2613722.0620459,50438504.9828027,577561.228435556,698939.061389366,120591.598619996,813961.887834364,1058150.45418467,11742861.0460984,4356370.91861309,10631749.1028242]) #Putman_2021 calculation
M31Sat_m=np.array([122462743.288989,1637600.50593105,33257712.3203104,4821095.889461,74198190.1960031,17360575.1453082,1029003.9565436,5471716.43137707,1496320.24410118,2673541.00325403,32526318.994349,19273771.7204746,45752990.31932,7638624.77810868,61011784.3039265,30787433.8275809,24132225.3714772,13646376.2387283,19380792.6353565,48189223.4697131,37483484.5985755,190992472.379195,5135817.58442046,433145333.818762,2968616.14441648,86753096.39298,336243115.244791,44124177.9592858,7941780.919495,2989197.80532057,919140.502147116,1018959.69662228])

M31Sat_sigma=np.array([24,4,7.8,6.4,16,11.8,2.9,5.7,4.6,7.1,10.2,9.3,8.6,10.9,7.1,11.5,8.4,5.3,92,3,7.8,35,14.8,4.5,5.8,2.6])

def PlotShpere(ax,R,xc,yc,zc):
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = xc+(R/1000.)*np.cos(u)*np.sin(v)
	y = yc+(R/1000.)*np.sin(u)*np.sin(v)
	z = zc+(R/1000.)*np.cos(v)
	ax.plot_wireframe(x, y, z, color="gray",alpha=0.3)


#how to run: python SelectHalo.py halo_catalog_1 halo_catalog_2 num_limit_for_dwarfs M_target_max M_target_min M_dwarf v_min
#example: $python distributions.py halos_0.0_G.ascii halos_0.0_C.ascii 250 1.3e12 1.1e12 1.0e8 4

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("HaloCatalogG",type=str)
    parser.add_argument("HaloCatalogC",type=str)
    parser.add_argument("NumLimit", type=int)
    parser.add_argument("M_highLim", type=float)
    parser.add_argument("M_lowLim", type=float)
    parser.add_argument("M_DwarfLim", type=float)
    parser.add_argument("V_min", type=float)
    args = parser.parse_args()
    dg=np.genfromtxt(args.HaloCatalogG, skip_header=18)#,names=True, skip_header=5)
    dc=np.genfromtxt(args.HaloCatalogC,skip_header=18)#names=True, skip_header=5)
    IDGAll=np.array(dg[:,0])
    IDCAll=np.array(dc[:,0])
    NumGAll=np.array(dg[:,1])
    NumCAll=np.array(dc[:,1])
    #print len(dc[:,2])-len(dg[:,2]) # mvir
    MgAll=np.array(dg[:,2])
    McAll=np.array(dc[:,2])
    XgAll=np.array(dg[:,8])
    YgAll=np.array(dg[:,9])
    ZgAll=np.array(dg[:,10])
    XcAll=np.array(dc[:,8])
    YcAll=np.array(dc[:,9])
    ZcAll=np.array(dc[:,10])
    RvGAll=np.array(dg[:,4])
    RvCAll=np.array(dc[:,4])
    VmaxGAll=np.array(dg[:,5])
    VmaxCAll=np.array(dc[:,5])
    ## Criterias
    MHaloHigh=args.M_highLim
    MHaloLow=args.M_lowLim
    ProbeRadius=0.3 # kpc -> Mpc
    NLim=args.NumLimit
    MDLim=args.M_DwarfLim
    V_min=args.V_min
    nbins=60
    nbins2=17
    ## Mass Sample
    MG=MgAll[NumGAll>NLim]
    MC=McAll[NumCAll>NLim]
    XG=XgAll[NumGAll>NLim]
    YG=YgAll[NumGAll>NLim]
    ZG=ZgAll[NumGAll>NLim]
    XC=XcAll[NumCAll>NLim]
    YC=YcAll[NumCAll>NLim]
    ZC=ZcAll[NumCAll>NLim]
    RvG=(RvGAll[NumGAll>NLim])/1000.
    RvC=(RvCAll[NumCAll>NLim])/1000.
    IDG=IDGAll[NumGAll>NLim]
    IDC=IDCAll[NumCAll>NLim]
    VmaxG=VmaxGAll[NumGAll>NLim]
    VmaxC=VmaxCAll[NumCAll>NLim]
    #Find target halos
    IDHaloG=IDG[(MG>MHaloLow) & (MG<MHaloHigh)]
    MHaloG=MG[(MG>MHaloLow) & (MG<MHaloHigh)]
    XHaloG=XG[(MG>MHaloLow) & (MG<MHaloHigh)]
    YHaloG=YG[(MG>MHaloLow) & (MG<MHaloHigh)]
    ZHaloG=ZG[(MG>MHaloLow) & (MG<MHaloHigh)]
    RvHaloG=RvG[(MG>MHaloLow) & (MG<MHaloHigh)]
    IDHaloC=IDC[(MC>MHaloLow) & (MC<MHaloHigh)]
    MHaloC=MC[(MC>MHaloLow) & (MC<MHaloHigh)]
    XHaloC=XC[(MC>MHaloLow) & (MC<MHaloHigh)]
    YHaloC=YC[(MC>MHaloLow) & (MC<MHaloHigh)]
    ZHaloC=ZC[(MC>MHaloLow) & (MC<MHaloHigh)]
    RvHaloC=RvC[(MC>MHaloLow) & (MC<MHaloHigh)]
    #print some info
    print("No of the main halos:")
    print("DMO:",len(MHaloG))
    print("CoSANG:",len(MHaloC))
    #let's find the most massive halo in the given range:
    MMG_index=np.argmax(MHaloG)
    MMC_index=np.argmax(MHaloC)
    #Now let's extract the dwarfs
    IDDwarfG=IDG[MG<MDLim]
    MDwarfG=MG[MG<MDLim]
    XDwarfG=XG[MG<MDLim]
    YDwarfG=YG[MG<MDLim]
    ZDwarfG=ZG[MG<MDLim]
    RvDwarfG=RvG[MG<MDLim]
    IDDwarfC=IDC[MC<MDLim]
    MDwarfC=MC[MC<MDLim]
    XDwarfC=XC[MC<MDLim]
    YDwarfC=YC[MC<MDLim]
    ZDwarfC=ZC[MC<MDLim]
    RvDwarfC=RvC[MC<MDLim]
    VmaxDwarfG=VmaxG[MG<MDLim]
    VmaxDwarfC=VmaxC[MC<MDLim]
    #################### v max limit
    IDDwarfG=IDDwarfG[VmaxDwarfG>V_min]
    MDwarfG=MDwarfG[VmaxDwarfG>V_min]
    XDwarfG=XDwarfG[VmaxDwarfG>V_min]
    YDwarfG=YDwarfG[VmaxDwarfG>V_min]
    ZDwarfG=ZDwarfG[VmaxDwarfG>V_min]
    RvDwarfG=RvDwarfG[VmaxDwarfG>V_min]
    IDDwarfC=IDDwarfC[VmaxDwarfC>V_min]
    MDwarfC=MDwarfC[VmaxDwarfC>V_min]
    XDwarfC=XDwarfC[VmaxDwarfC>V_min]
    YDwarfC=YDwarfC[VmaxDwarfC>V_min]
    ZDwarfC=ZDwarfC[VmaxDwarfC>V_min]
    RvDwarfC=RvDwarfC[VmaxDwarfC>V_min]#
    VmaxDwarfG=VmaxDwarfG[VmaxDwarfG>V_min]
    VmaxDwarfC=VmaxDwarfC[VmaxDwarfC>V_min]
	####################
    #print len(MDwarfG)
    #print("No of halos belove %.4g"%DHL,"is:")
    #print("DMO:",len(MDwarfG))
    #print("NoDisk:",len(MDwarfC0))
    #print("CoSANG:",len(MDwarfC))
    ###############
    # Let's do analysis for the most massive halo first
    XHG=XHaloG[MMG_index]
    YHG=YHaloG[MMG_index]
    ZHG=ZHaloG[MMG_index]
    RvHG=RvHaloG[MMG_index]
    XHC=XHaloC[MMC_index]
    YHC=YHaloC[MMC_index]
    ZHC=ZHaloC[MMC_index]
    RvHC=RvHaloC[MMC_index]
    print("Halo mass:%g"%MHaloG[MMG_index])
    print("Halo Rv:%g"%RvHaloG[MMG_index])
    print("Halo ID:%g"%IDHaloC[MMC_index])
    print("Halo pos:%g,%g,%g"%(XHaloC,YHaloC,ZHaloC))
    #let's plot these
    rHDwarfsG=np.sqrt((XDwarfG-XHG)**2.0+(YDwarfG-YHG)**2.0+(ZDwarfG-ZHG)**2.0)
    rHDwarfsC=np.sqrt((XDwarfC-XHC)**2.0+(YDwarfC-YHC)**2.0+(ZDwarfC-ZHC)**2.0)
    rGlim=1.0*RvHG
    rClim=1.0*RvHC
    xg1=XDwarfG[rHDwarfsG<rGlim]
    yg1=YDwarfG[rHDwarfsG<rGlim]
    zg1=ZDwarfG[rHDwarfsG<rGlim]
    mg1=MDwarfG[rHDwarfsG<rGlim]
    xc1=XDwarfC[rHDwarfsC<rClim]
    yc1=YDwarfC[rHDwarfsC<rClim]
    zc1=ZDwarfC[rHDwarfsC<rClim]
    mc1=MDwarfC[rHDwarfsC<rClim]
    rg1=rHDwarfsG[rHDwarfsG<rGlim]
    rc1=rHDwarfsC[rHDwarfsC<rClim]
    #print some info
    print("No of the sub halos:")
    print("DMO:",len(mg1))
    print("CoSANG:",len(mc1))
    #now let's extract all info for all halos in the given range
    #
    rG=[]
    #rg=[None]*len(IDHaloG)
    for idg in IDHaloG:
        rg2=np.sqrt((XHaloG[IDHaloG==idg]-XDwarfG)**2.+(YHaloG[IDHaloG==idg]-YDwarfG)**2.+(ZHaloG[IDHaloG==idg]-ZDwarfG)**2. )
        rG.extend(rg2)#[rg2<(RvHaloG[IDHaloG==idg]/1000.)])
    rC=[]
    #rg=[None]*len(IDHaloG)
    #for idc in IDHaloC:
    #    rc2=np.sqrt((XHaloC[IDHaloC==idc]-XDwarfC)**2.+(YHaloC[IDHaloC==idc]-YDwarfC)**2.+(ZHaloC[IDHaloC==idc]-ZDwarfC)**2. )
    #    rC.extend(rc2)#[rc2<(RvHaloC[IDHaloC==idc]/1000.)])
    #
    #rGArray=np.array(rG)
    #rCArray=np.array(rC)
    #print(rG)
    #rG2=rG#[rGArray<0.2]
    #rC2=rC#Array[rCArray<0.2]
    #plots
    fig = plt.figure(1)
    #fig.suptitle('CoSANG vs N-Body ')
    ax11 = fig.add_subplot(221)
    ax11.set_xlabel('X (Mpc)')
    ax11.set_ylabel('Y (Mpc)')
    ax11.set_title('Subhalo Distribution G')
    #ax11.plot(XHG,YHG,'b+')
    ax11.scatter(XHG,YHG,c='red', alpha=0.99, marker='o',s=60)#RvHG)
    ax11.scatter(xg1,yg1,c='black', alpha=0.7, marker='.',s=15)
    ax12 = fig.add_subplot(222)
    ax12.set_xlabel('X (Mpc)')
    ax12.set_ylabel('Y (Mpc)')
    ax12.set_title('Subhalo Distribution C')
    ax12.scatter(XHC,YHC,c='red', alpha=0.99, marker='o',s=60)#RvHC)
    ax12.scatter(xc1,yc1,c='black', alpha=0.7, marker='.',s=15)
    ax13 = fig.add_subplot(223)
    ax13.set_xlabel('X (Mpc)')
    ax13.set_ylabel('Z (Mpc)')
    ax13.set_title('Subhalo Distribution G')
    ax13.scatter(XHG,ZHG,c='red', alpha=0.99, marker='o',s=60)#RvHG)
    ax13.scatter(xg1,zg1,c='black', alpha=0.7, marker='.',s=15)
    ax14 = fig.add_subplot(224)
    ax14.set_xlabel('X (Mpc)')
    ax14.set_ylabel('Z (Mpc)')
    ax14.set_title('Subhalo Distribution C')
    ax14.scatter(XHC,ZHC,c='red', alpha=0.99, marker='o',s=60)#RvHC)
    ax14.scatter(xc1,zc1,c='black', alpha=0.7, marker='.',s=15)
    #histogram of aboundances
    fig2 = plt.figure(2,figsize=plt.figaspect(2))
    #fig2.suptitle('CoSANG vs N-Body ')
    ax21 = fig2.add_subplot(211)
    ax21.set_xlabel('$Log (M_{halo})[M_{\odot}]$')
    ax21.set_ylabel('$N_M(>M)$')
    ax21.set_title('Dwarf Halos Mass Abundance')
    ax21.set_xlim(6.0,10.8)
    #ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
    ax21.hist(np.log10(mg1),linewidth=2, bins=nbins, log=True,cumulative=-1, histtype='step', alpha=0.9,color='blue',label='DMO')
    ax21.hist(np.log10(mc1),linewidth=2,bins=nbins,log=True, cumulative=-1, histtype='step', alpha=0.9,color='red',label='CoSANG')
    ax21.hist(np.log10(MWSat_m),linewidth=2,linestyle='-.',bins=nbins2,log=True, cumulative=-1, histtype='step', alpha=0.9,color='black',label='MW')
    ax21.hist(np.log10(M31Sat_m),linewidth=2,linestyle='-.',bins=nbins2,log=True, cumulative=-1, histtype='step', alpha=0.9,color='gray',label='M31')
    ax21.legend(loc=1)
    ax22 = fig2.add_subplot(212)
    ax22.set_xlabel('$r_{halo}/R_v$')
    ax22.set_ylabel('$N_r(<d)$')
    ax22.set_title('Dwarf Halos Distance Abundance')
    ax22.set_xlim(0.06,0.98)
    #ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
    ax22.hist(rg1/RvHG,linewidth=2, bins=nbins, log=True,cumulative=True, histtype='step', alpha=0.9,color='blue',label='DMO')
    ax22.hist(rc1/RvHC,linewidth=2,bins=nbins,log=True,cumulative=True, histtype='step', alpha=0.9,color='red',label='CoSANG')
    MWSat2=MWSat/140
    ax22.hist(MWSat2,linewidth=2,linestyle='-.',bins=nbins,log=True,cumulative=True, histtype='step', alpha=0.9,color='black',label='MW')
    M31Sat2=M31Sat/140
    ax22.hist(M31Sat2,linewidth=2,linestyle='-.',bins=nbins,log=True,cumulative=True, histtype='step', alpha=0.9,color='gray',label='M31')
    ax22.legend(loc=2)
    #fig3 = plt.figure(3)
    #ax3 = fig3.add_subplot(111)
    #ax3.set_xlabel('$r_{halo}[kpc]$')
    #ax3.set_ylabel('$N$')
    #ax3.set_title('Total Dwarf Halos Distance aboundance')
    #ax3.hist(rG2*1000,linewidth=2, bins=nbins2, log=False, histtype='step', alpha=0.9,color='blue',label='DMO')
    #ax3.hist(rC2*1000,linewidth=2,bins=nbins2,log=False, histtype='step', alpha=0.9,color='green',label='CoSANG')
    #ax3.legend(loc=2)
    #We are done with the most massive halo, now let's look at all massive halos
    fig3 = plt.figure(3,figsize=plt.figaspect(0.5))
    ax31=fig3.add_subplot(121) #log log + cumulative top down,
    ax31.set_xlabel('$Log(V_{max})$')
    ax31.set_ylabel('$N_{V_{max}}$')
    ax31.title.set_text(h_title)
    ax31.hist(np.log10(VmaxDwarfG),linewidth=2, bins=nbins, log=True,cumulative=False, histtype='step', alpha=0.9,color='blue',label='DMO')
    ax31.hist(np.log10(VmaxDwarfC),linewidth=2, bins=nbins, log=True,cumulative=False, histtype='step', alpha=0.9,color='red',label='CoSANG')
    ax31.hist(np.log10(MWSat_sigma),linewidth=2,linestyle='-.', bins=nbins2, log=True,cumulative=False, histtype='step', alpha=0.9,color='black',label='MW')
    ax31.hist(np.log10(M31Sat_sigma),linewidth=2,linestyle='-.', bins=nbins2, log=True,cumulative=False, histtype='step', alpha=0.9,color='gray',label='M31')
    ax31.set_xlim(0.06,1.9)
    ax31.set_ylim(0.7,1500)
    ax31.legend(loc=1)
    ax32=fig3.add_subplot(122)
    ax32.set_xlabel('$Log(V_{max})$')
    ax32.set_ylabel('$N(>V_{max})$') # np.log10(VmaxDwarfG)
    ax32.title.set_text(h_title)
    ax32.set_xlim(0.06,1.9)
    ax32.set_ylim(0.7,10000)
    ax32.hist(np.log10(VmaxDwarfG),linewidth=2, bins=nbins, log=True,cumulative=-1, histtype='step', alpha=0.9,color='blue',label='DMO')
    ax32.hist(np.log10(VmaxDwarfC),linewidth=2, bins=nbins, log=True,cumulative=-1, histtype='step', alpha=0.9,color='red',label='CoSANG')
    ax32.hist(np.log10(MWSat_sigma),linewidth=2,linestyle='-.', bins=nbins2, log=True,cumulative=-1, histtype='step', alpha=0.9,color='black',label='MW')
    ax32.hist(np.log10(M31Sat_sigma),linewidth=2,linestyle='-.', bins=nbins2, log=True,cumulative=-1, histtype='step', alpha=0.9,color='gray',label='M31')
    Kelly_v=[0.653,1.54]
    Kelly_N=[2000,1]
    ax32.plot(Kelly_v,Kelly_N,linewidth=2,linestyle=':',label='Kelly+19',color='gray')
    ax32.legend(loc=1)
    plt.show()


'''

#### Find most massive halo
MMG_index=np.argmin(MMassiveG)
MMC_index=np.argmin(MMassiveC)
MMC0_index=np.argmin(MMassiveC0)

rMg=np.sqrt(XMassiveG[MMG_index]**2. + YMassiveG[MMG_index]**2. + ZMassiveG[MMG_index]**2.)
rMc0=np.sqrt(XMassiveC0[MMC0_index]**2. + YMassiveC0[MMC0_index]**2. + ZMassiveC0[MMC0_index]**2.)
rMc=np.sqrt(XMassiveC[MMC_index]**2. + YMassiveC[MMC_index]**2. + ZMassiveC[MMC_index]**2.)


print("We selected this massive halo:")
print("DMO=%.4g"%MMassiveG[MMG_index],"solar mass & coordinates:",XMassiveG[MMG_index],",",YMassiveG[MMG_index],",",ZMassiveG[MMG_index])
print("NoDisk=%.4g"%MMassiveC0[MMC0_index],"solar mass & coordinates:",XMassiveC0[MMC0_index],",",YMassiveC0[MMC0_index],",",ZMassiveC0[MMC0_index])
print("CoSANG=%.4g"%MMassiveC[MMC_index],"solar mass & coordinates:",XMassiveC[MMC_index],",",YMassiveC[MMC_index],",",ZMassiveC[MMC_index])


print("No of sattelites closer than",ProbeRadius*1000," kpc is:")
print("DMO:",len(rSampleG))
print("NoDisk:",len(rSampleC0))
print("CoSANG:",len(rSampleC))



########################################################## plots
fig = plt.figure(1)
fig.suptitle('CoSANG vs N-Body ')

ax1 = fig.add_subplot(221)
ax1.set_xlabel('$Log (M_{halo})$')
ax1.set_ylabel('$N$')
ax1.set_title('Dwarf Halos Mass aboundance')
#ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
ax1.hist(np.log10(MDwarfG),linewidth=2, bins=nbins, log=False, histtype='step', alpha=0.9,color='blue',label='DMO')
ax1.hist(np.log10(MDwarfC0),linewidth=2,bins=nbins,log=False, histtype='step', alpha=0.9,color='black',label='NoDisk')
ax1.hist(np.log10(MDwarfC),linewidth=2,bins=nbins,log=False, histtype='step', alpha=0.9,color='green',label='CoSANG')

ax1.legend(loc=2)


#########################################################
ax2 = fig.add_subplot(222)

ax2.set_xlabel('$Log (M_{halo})$')
ax2.set_ylabel('$N$')
ax2.set_title('Mass aboundance')
colors=['green','black','blue']
labels=['CoSANG','NoDisk','DMO']

ax2.hist([Mclog, Mc0log, Mglog], bins=nbins, log=True, color=colors,label=labels)
ax2.legend()
#########################################################
ax3 = fig.add_subplot(223)
ax3.set_xlabel('$log(d_{kpc})$')
ax3.set_ylabel('$N$')
ax3.set_title('Subhalo Distribution for ~%.1g halo'%MMassiveG[MMG_index])

#Mc=np.log10(dc[:,2])
#Mg=np.log10(dg[:,2])



ng, binsg, patchesg= ax3.hist(np.log10(rSampleG),linewidth=2,bins=nbins,log=False,cumulative=True, histtype='step', alpha=0.8,color='blue',label='DMO')
nc0, binsc0, patchesc0= ax3.hist(np.log10(rSampleC0),linewidth=2,bins=nbins,log=False,cumulative=True,histtype='step', alpha=0.8,color='black',label='NoDisk')
nc, binsc, patchesc= ax3.hist(np.log10(rSampleC),linewidth=2,bins=nbins,log=False, cumulative=True,histtype='step', alpha=0.8,color='green',label='CoSANG')
#ax3.legend()


#ng, binsg, patchesg= ax3.hist(rSampleG,bins=100,log=True, alpha=0.6,color='green',label='NoDisk')
#nc, binsc, patchesc= ax3.hist(rSampleC,bins=100,log=True, alpha=0.4,color='blue',label='CoSANG')
#ax3.legend()
ax3.legend( loc =2)


#sns.distplot(np.log10(rSampleC), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'CoSANG')
#sns.distplot(np.log10(rSampleG), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'NoDisk')

#######################################################

ax4=fig.add_subplot(224)
ax4.set_xlabel('$Bin number~d_{kpc}$')
ax4.set_ylabel('$n_i/n_g$')
ax4.set_title('Subhalo Distribution')
colors=['grey','black']
labels=['CoSANG','NoDisk']


#######################################################
fig2 = plt.figure(2)
#fig2.set_title('Mass distribution')
#sns.kdeplot(Mc)
ax2_1=fig2.add_subplot(111)
ax2_1.set_xlabel('$Log(M)$')
ax2_1.set_ylabel('$n$')
ax2_1.set_title('Dwarf Halos Mass Distribution')
sns.distplot(np.log10(MDwarfC), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'CoSANG')
sns.distplot(np.log10(MDwarfC0), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'NoDisk')
sns.distplot(np.log10(MDwarfG), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'DMO')
#fig2.legend(prop={'size': 16}, title = 'Airline')

###########################
plt.show()
'''
