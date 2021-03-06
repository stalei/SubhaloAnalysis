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

def PlotShpere(ax,R,xc,yc,zc):
	u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
	x = xc+(R/1000.)*np.cos(u)*np.sin(v)
	y = yc+(R/1000.)*np.sin(u)*np.sin(v)
	z = zc+(R/1000.)*np.cos(v)
	ax.plot_wireframe(x, y, z, color="gray",alpha=0.3)


#how to run: python SelectHalo.py halo_catalog_1 halo_catalog_2 num_limit_for_dwarfs M_target_max M_target_min M_dwarf
#example: $python distributions.py halos_0.0_G.ascii halos_0.0_C.ascii 250 1.3e12 1.1e12 1.0e8

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("HaloCatalogG",type=str)
    parser.add_argument("HaloCatalogC",type=str)
    parser.add_argument("NumLimit", type=int)
    parser.add_argument("M_highLim", type=float)
    parser.add_argument("M_lowLim", type=float)
    parser.add_argument("M_DwarfLim", type=float)
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
    ## Criterias
    MHaloHigh=args.M_highLim
    MHaloLow=args.M_lowLim
    ProbeRadius=0.3 # kpc -> Mpc
    NLim=args.NumLimit
    MDLim=args.M_DwarfLim
    nbins=30
    nbins2=60
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
    print("Gadget:",len(MHaloG))
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
    #print len(MDwarfG)
    #print("No of halos belove %.4g"%DHL,"is:")
    #print("Gadget:",len(MDwarfG))
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
    print("Gadget:",len(mg1))
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
    ax11.set_xlabel('X')
    ax11.set_ylabel('Y')
    ax11.set_title('Dwarfs Distribution G')
    #ax11.plot(XHG,YHG,'b+')
    ax11.scatter(XHG,YHG,c='black', alpha=0.9, marker='+',s=50)#RvHG)
    ax11.scatter(xg1,yg1,c='black', alpha=0.7, marker='o',s=15)
    ax12 = fig.add_subplot(222)
    ax12.set_xlabel('X')
    ax12.set_ylabel('Y')
    ax12.set_title('Dwarfs Distribution C')
    ax12.scatter(XHC,YHC,c='black', alpha=0.9, marker='+',s=50)#RvHC)
    ax12.scatter(xc1,yc1,c='black', alpha=0.7, marker='o',s=15)
    ax13 = fig.add_subplot(223)
    ax13.set_xlabel('X')
    ax13.set_ylabel('Z')
    ax13.set_title('Dwarfs Distribution G')
    ax13.scatter(XHG,ZHG,c='black', alpha=0.9, marker='+',s=50)#RvHG)
    ax13.scatter(xg1,zg1,c='black', alpha=0.7, marker='o',s=15)
    ax14 = fig.add_subplot(224)
    ax14.set_xlabel('X')
    ax14.set_ylabel('Z')
    ax14.set_title('Dwarfs Distribution C')
    ax14.scatter(XHC,ZHC,c='black', alpha=0.9, marker='+',s=50)#RvHC)
    ax14.scatter(xc1,zc1,c='black', alpha=0.7, marker='o',s=15)
    #histogram of aboundances
    fig2 = plt.figure(2,figsize=plt.figaspect(2))
    #fig2.suptitle('CoSANG vs N-Body ')
    ax21 = fig2.add_subplot(211)
    ax21.set_xlabel('$Log (M_{halo})$')
    ax21.set_ylabel('$N_M(>M)$')
    ax21.set_title('Dwarf Halos Mass abundance')
    #ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
    ax21.hist(np.log10(mg1),linewidth=2, bins=nbins, log=False,cumulative=-1, histtype='step', alpha=0.9,color='blue',label='Gadget')
    ax21.hist(np.log10(mc1),linewidth=2,bins=nbins,log=False, cumulative=-1, histtype='step', alpha=0.9,color='red',label='CoSANG')
    ax21.legend(loc=1)
    ax22 = fig2.add_subplot(212)
    ax22.set_xlabel('$r_{halo}[kpc]$')
    ax22.set_ylabel('$N_r(<d)$')
    ax22.set_title('Dwarf Halos Distance abundance')
    #ax1.plot(d1['star_age'],d1['center h1']) #plot of main data
    ax22.hist(rg1*1000,linewidth=2, bins=nbins, log=False,cumulative=True, histtype='step', alpha=0.9,color='blue',label='Gadget')
    ax22.hist(rc1*1000,linewidth=2,bins=nbins,log=False,cumulative=True, histtype='step', alpha=0.9,color='red',label='CoSANG')
    ax22.legend(loc=2)
    #fig3 = plt.figure(3)
    #ax3 = fig3.add_subplot(111)
    #ax3.set_xlabel('$r_{halo}[kpc]$')
    #ax3.set_ylabel('$N$')
    #ax3.set_title('Total Dwarf Halos Distance aboundance')
    #ax3.hist(rG2*1000,linewidth=2, bins=nbins2, log=False, histtype='step', alpha=0.9,color='blue',label='Gadget')
    #ax3.hist(rC2*1000,linewidth=2,bins=nbins2,log=False, histtype='step', alpha=0.9,color='green',label='CoSANG')
    #ax3.legend(loc=2)
    #We are done with the most massive halo, now let's look at all massive halos
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
print("Gadget=%.4g"%MMassiveG[MMG_index],"solar mass & coordinates:",XMassiveG[MMG_index],",",YMassiveG[MMG_index],",",ZMassiveG[MMG_index])
print("NoDisk=%.4g"%MMassiveC0[MMC0_index],"solar mass & coordinates:",XMassiveC0[MMC0_index],",",YMassiveC0[MMC0_index],",",ZMassiveC0[MMC0_index])
print("CoSANG=%.4g"%MMassiveC[MMC_index],"solar mass & coordinates:",XMassiveC[MMC_index],",",YMassiveC[MMC_index],",",ZMassiveC[MMC_index])


print("No of sattelites closer than",ProbeRadius*1000," kpc is:")
print("Gadget:",len(rSampleG))
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
ax1.hist(np.log10(MDwarfG),linewidth=2, bins=nbins, log=False, histtype='step', alpha=0.9,color='blue',label='Gadget')
ax1.hist(np.log10(MDwarfC0),linewidth=2,bins=nbins,log=False, histtype='step', alpha=0.9,color='black',label='NoDisk')
ax1.hist(np.log10(MDwarfC),linewidth=2,bins=nbins,log=False, histtype='step', alpha=0.9,color='green',label='CoSANG')

ax1.legend(loc=2)


#########################################################
ax2 = fig.add_subplot(222)

ax2.set_xlabel('$Log (M_{halo})$')
ax2.set_ylabel('$N$')
ax2.set_title('Mass aboundance')
colors=['green','black','blue']
labels=['CoSANG','NoDisk','Gadget']

ax2.hist([Mclog, Mc0log, Mglog], bins=nbins, log=True, color=colors,label=labels)
ax2.legend()
#########################################################
ax3 = fig.add_subplot(223)
ax3.set_xlabel('$log(d_{kpc})$')
ax3.set_ylabel('$N$')
ax3.set_title('Subhalo Distribution for ~%.1g halo'%MMassiveG[MMG_index])

#Mc=np.log10(dc[:,2])
#Mg=np.log10(dg[:,2])



ng, binsg, patchesg= ax3.hist(np.log10(rSampleG),linewidth=2,bins=nbins,log=False,cumulative=True, histtype='step', alpha=0.8,color='blue',label='Gadget')
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
sns.distplot(np.log10(MDwarfG), hist = False, kde = True,kde_kws = {'linewidth': 3},label = 'Gadget')
#fig2.legend(prop={'size': 16}, title = 'Airline')

###########################
plt.show()
'''
