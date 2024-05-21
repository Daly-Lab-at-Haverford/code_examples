#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import psi4
import v3d


# In[ ]:


psi4.core.set_output_file('bonds.dat')
theory = 'mp2/6-31G*'
bohr2ang = 0.529177
rad2deg = 180.0/np.pi
#psi4.set_memory('2gb')


# In[ ]:


#optimized coordinates of the molecules
azi = """
0 1
    C           -1.372186183321    -2.010860020902    -0.886820192828
    C           -0.173158673000    -1.298344357907    -0.406115006426
    N            0.786309071860    -1.384810968974     0.420697487368
    C            0.662837700048    -0.111043752576    -0.491918176988
    H            1.422532833369    -0.002132476819    -1.261568419806
    C            0.271415655367     1.154529325936     0.168585840901
    O           -0.638034693748     0.950769479245     1.162558403008
    O            0.686902544636     2.258017066085    -0.144478115352
    H           -1.457401352716    -2.978238526949    -0.387585199421
    H           -1.313874303487    -2.157164437091    -1.968871402470
    H           -2.260422649854    -1.409349289546    -0.673788389312
    H           -0.815494634684     1.839688017493     1.533442847642
"""


# In[ ]:


#what coordinates are we optimizing
coords = [
    (2, 4),
    (6, 4),
    (4, 3)
]

opt_coords = np.zeros(len(coords))
#print the optimized values of these coordinates:
#also, save these values for later
azi_geo = psi4.geometry(azi)
for item, pair in enumerate(coords):
    opt_coords[item] = v3d.dist((azi_geo.to_arrays()[0])[pair[0]-1],
                 (azi_geo.to_arrays()[0])[pair[1]-1])*bohr2ang
    print(opt_coords[item])

opt_energy = psi4.energy(theory)


# In[ ]:


#collecting real optimized coordinate values here. 
#these wont always be the same as the values in r
#because of "fixed" rather than "frozen" coordinate scan
r_real = np.zeros((len(coords), 15))
#collecting energies here
energies = np.zeros((len(coords), 15))
    
for item, pair in enumerate(coords):
    
    #define geometry
    azi_geo = psi4.geometry(azi)
    #get optimized coordinate values
    r0 = v3d.dist((azi_geo.to_arrays()[0])[pair[0]-1],
                  (azi_geo.to_arrays()[0])[pair[1]-1])*bohr2ang
    
    #increased coordinate value half of scan
    rplus = np.linspace(r0, r0+0.15, 8)
    #decreased coordinate value half of scan
    rminus = np.linspace(r0-0.15, r0, 8)[:-1]
    #putting them together
    r = np.append(rminus, rplus)
        
    azi_geo = psi4.geometry(azi)
    #selecting left half. This set up 
    #lets the optimization start with the optimized
    #coordinates, then move away from them
    #hence, np.flip.
    for i, n in enumerate(np.flip(rminus)):#r):
        #options for geometry opt
        options={'scf_type':'df',
             'g_convergence':'gau_loose',
             'freeze_core':'true',
             'opt_coordinates':'both',
             #selecting atoms I care about and target distance
             "ranged_distance":"{a} {b} {c} {d}".format(a=pair[0], 
                                                   b=pair[1], 
                                                   c=n-0.005,
                                                   d=n+0.005),
             "geom_maxiter" : 200,
             #"full_hess_every" : 20,
             "dynamic_level" : 0}
        psi4.set_options(options)
        energies[item, 6-i] = psi4.optimize(theory)
        r_real[item, 6-i] = v3d.dist((azi_geo.to_arrays()[0])[pair[0]-1],
                               (azi_geo.to_arrays()[0])[pair[1]-1])*bohr2ang
        print(np.abs(r_real[item, 6-i] - n))
            
    azi_geo = psi4.geometry(azi)
    #generally the same as above but for the above half of the scan
    for i, n in enumerate(rplus):#r):
        options={'scf_type':'df',
             'g_convergence':'gau_loose',
             'freeze_core':'true',
             'opt_coordinates':'both',
             "ranged_distance":"{a} {b} {c} {d}".format(a=pair[0], 
                                                   b=pair[1], 
                                                   c=n-0.005,
                                                   d=n+0.005),             
             "geom_maxiter" : 200,
             #"full_hess_every" : 20,
             "dynamic_level" : 0}
        psi4.set_options(options)
        energies[item, 7+i] = psi4.optimize(theory)
        r_real[item, 7+i] = v3d.dist((azi_geo.to_arrays()[0])[pair[0]-1],
                               (azi_geo.to_arrays()[0])[pair[1]-1])*bohr2ang
        print(np.abs(r_real[item, 7+i] - n))


# In[6]:


print(energies)


# In[14]:


energies_norm = (energies-opt_energy)*627.503


# In[17]:


#collected data plots
colors = ['r', 'g', 'b']

for i in range(len(coords)):
               plt.plot(r_real[i], energies_norm[i], '--o', color=colors[i]) 
               #plt.ylim(0, 10)
               #plt.xlim(1.0, 2)
               plt.show()

# In[18]:


def vbond(r, k, r0):
    return k*(r-r0)**2


# In[19]:


parms = []
                 
for i in range(len(coords)):
    parms.append(spo.curve_fit(vbond, 
                               r_real[i], 
                               energies_norm[i], 
                               p0=[300, 1.5]))


# In[20]:


#fitted data plots
colors = ['r', 'g', 'b']

for i in range(len(coords)):
    plt.plot(r_real[i], energies_norm[i], 'o', color=colors[i], label=str(coords[i])) 
    plt.plot(r_real[i], vbond(r_real[i], *parms[i][0]), '--', color=colors[i])
plt.legend()
plt.show()
#plt.ylim(0, 50)
#plt.xlim(1.0, 2)


# In[21]:


for i in range(len(coords)):
    print('For the bond between atoms ', str(coords[i]))
    print('the force constant is ', str(parms[i][0][0]), ' p/m ', np.sqrt(np.diag(parms[i][1]))[0], 'kcal/mol/A^2')
    print('and the optimal bond length is ', str(parms[i][0][1]), ' p/m ', np.sqrt(np.diag(parms[i][1]))[1], 'A')


# In[ ]:





# In[ ]:




