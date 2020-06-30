import numpy as np
import networkx as nx

from scipy.signal import savgol_filter
from scipy.signal import find_peaks

sizes = [100, 100, 1000,1000]
probs = [[0.3, 0.05, 0.2, 0.01],[0.07, 0.3, 0.3, 0.01],[0.01, 0.01, 0.1, 0.01],[0.005,0.005,0.03,0.01]]
p=np.ones(np.sum(sizes))
#p[np."random.ranint(len(p))]="-1 
p[3]=-1





def thread_func(tq,tinc=[5,10,15],n_run=10,n_spread=100):
    tot_peaks=np.zeros(len(tinc))
    for k in range(0,len(tinc)):
        print(tq,tinc[k],"START")
        for j in range(0,n_run):
            covid=pandemic(sizes,probs,p,incubation_time=tinc[k],qtime=tq)
            for i in range(0,n_spread): #spread 100 times
                covid.spread()
                
            S,E,R,C,D=covid.get_SIR_data() #final data
            xf=savgol_filter(E,5,2) #filter
            peaks, _ = find_peaks(xf, distance=10) #find peaks
            tot_peaks[k]+=len(peaks)
        print(tq,tinc[k],tot_peaks)
        tot_peaks[k]/=n_run
        print(tq,tinc[k],"STOP")
        
    return tot_peaks

def thread_func2(titq,n_run=10,n_spread=100):
    tot_peaks=0
    print(titq,"START")
    for j in range(0,n_run):
        covid=pandemic(sizes,probs,p,incubation_time=titq[0],qtime=titq[1])
        for i in range(0,n_spread): #spread 100 times
            covid.spread()
                
        S,E,R,C,D=covid.get_SIR_data() #final data
        xf=savgol_filter(E,5,2) #filter
        peaks, _ = find_peaks(xf, distance=10) #find peaks
        tot_peaks+=len(peaks)
    tot_peaks/=n_run
    print(titq,"STOP")
        
    return tot_peaks

"""Pandemic class below..."""
class pandemic:
    def __init__(self,sizes,density,p0,directed=True,
                 beta=0.5,gamma=0.1,alpha=0.4,pi=0.5,mu=0.05,incubation_time=10,
                 ksi=0.01,delta=0.01,rho=0.005,nu=0.06,
                 qtime=20,use_quarantine=True,test_precision=1,social_tracking=1):
        """
        LOLOLOLLO
        adsa
        asdsad
        """
        self.sizes=sizes #size of blocks [N1,N2,N3,..] (elements qre no of nodes)
        self.density=density #density matrix (density of connections between blocks)
        self.plague=np.array([p0]) #stores the plague value of each node (1:S,-1:E,0:R,-2:C,-3:D)
        self.incubation_time=incubation_time
        
        self.quarantine=np.zeros((1,len(p0))) #quarantine timer
        self.incubation=np.zeros((1,len(p0))) #incubation timer
        self.incubation[0][np.ravel(np.argwhere(self.plague==-1))]=self.incubation_time
        self.use_quarantine=use_quarantine #if we want to implement quarantine
        self.qtime=qtime #lenght of quarantine
        
        self.beta=beta #p(E or C infects one of its neighbor)
        
        #rates for E
        #the next three HAVE TO SUM UP TO 1
        self.gamma=gamma #p(E->R|t>t_inc)
        self.alpha=alpha #p(E->C|t>t_inc)
        self.pi=pi #p(E->S|t>t_inc)
        self.mu=mu #probability that an exposed person gets tested for quarantine
        
        #rates for C
        self.ksi=ksi #p(C->D)
        self.delta=delta #p(C->R)
        self.rho=rho #p(C->S)
        self.nu=nu #probability that a critical person gets tested for quarantine
        
        self.test_precision=test_precision #p(positive test|infected)
        
        self.social_tracking=social_tracking
        #if we find an infected patient this is the probability that we find connection
        #with its i-th neighbouring node
        
        #to put only tested node to quarantine set this to 0
        
        #create graph
        self.G=nx.stochastic_block_model(sizes,probs,directed=directed,seed=123)
        self.age=0 #age of system
        
    
    def spread(self): #make one step infect or heal
        self.age+=1
        S,E,R,C,D=self.sort_ppl()
        self.plague=np.vstack([self.plague,self.plague[-1]])
        
        self.quarantine=np.vstack([self.quarantine,self.quarantine[-1]])
        self.quarantine[-1][np.ravel(np.argwhere(self.quarantine[-1]>0))]-=1 # step time in quarantine
        
        self.incubation=np.vstack([self.incubation,self.incubation[-1]])
        self.incubation[-1][np.ravel(np.argwhere(self.incubation[-1]>0))]-=1 # step incubation time
        
        #the last row will be the data we update
        #print(self.age)
        #handle exposed
        for i in range(0,len(E)):
            if self.quarantine[-1][i]==0:
                for j in [n for n in self.G.neighbors(E[i])]:
                    r=np.random.rand()
                    if self.plague[-1][j]==1 and r<self.beta and self.quarantine[-1][j]==0:
                        #with probability beta we infect each non-quarantined neighbor
                        self.plague[-1][j]=-1
                        #add random incubation period to newly infected node
                        self.incubation[-1][j]=np.round(np.random.exponential(self.incubation_time))
                if self.incubation[-1][i]==0:
                    #after incubation period patient can get resistant,susceptible or critical
                    state=np.random.choice([1,0,-2],p=[self.pi,self.gamma,self.alpha])
                    self.plague[-1][E[i]]=state
                    
                #quarantine
                #with probability mu we test patient with a "self.test_precise" precise test
                #if produces positive test we apply quarantine to the node and its neighbors
                if self.use_quarantine and self.plague[-1][E[i]]==-1:
                    r=np.random.rand()
                    if r<self.mu*self.test_precision: #with given probability and precision we test patient
                        self.quarantine[-1][i]=self.qtime #put given cell in quarainte
                        #try to track its neighbors aswell (with given success)
                        for j in [n for n in self.G.neighbors(E[i])]: 
                            if self.quarantine[-1][j]==0:
                                r=np.random.rand()
                                if r<self.social_tracking:
                                    self.quarantine[-1][j]=self.qtime
                                    
            else: #in quaraintine incubation time still steps
                if self.incubation[-1][i]==0:
                    #after incubation period patient can get resistant,susceptible or critical
                    state=np.random.choice([1,0,-2],p=[self.pi,self.gamma,self.alpha])
                    self.plague[-1][E[i]]=state
                    
            #handle critical
            for i in range(0,len(C)):
                if self.quarantine[-1][i]==0: #if not in quarantine they can still infect
                    #print("IN")
                    for j in [n for n in self.G.neighbors(C[i])]:
                        #print("IN j",j)
                        r=np.random.rand()
                        #print(self.plague[-1][j]==1,r<self.beta,self.quarantine[-1][j]==0)
                        if self.plague[-1][j]==1 and r<self.beta and self.quarantine[-1][j]==0:
                            #print("infect")
                            #with probability beta we infect each non-quarantined neighbor
                            self.plague[-1][j]=-1
                            #add random incubation period to newly infected node
                            self.incubation[-1][j]=np.round(np.random.exponential(self.incubation_time))
                
                #either in or not in quarantine C->R,D,S is possible
                r=np.random.rand(3)
                if r[0]<self.delta: #C->R
                    self.plague[-1][C[i]]=0
                elif r[1]<self.ksi: #C->D
                    self.plague[-1][C[i]]=-3
                elif r[2]<self.rho: #C->S
                    self.plague[-1][C[i]]=1
                
                        
                    
    
    def sort_ppl(self,pi=-1): #sort people to S,I,R
        S=np.ravel(np.argwhere(self.plague[pi]==1))
        E=np.ravel(np.argwhere(self.plague[pi]==-1))
        R=np.ravel(np.argwhere(self.plague[pi]==0))
        C=np.ravel(np.argwhere(self.plague[pi]==-2))
        D=np.ravel(np.argwhere(self.plague[pi]==-3))
        return S,E,R,C,D
    
    def reset(self): #reset the state
        self.age=0
        self.plague=np.array([self.plague[0]])
        self.quarantine=np.zeros(len(self.plague))
        self.incubation=np.zeros((1,len(self.plague))) #incubation timer
        self.incubation[0][np.ravel(np.argwhere(self.plague==-1))]=self.incubation_time

        
    def get_colors_shapes(self,p,q): #different color for each disease state, diff state for quarantine
        #, p is one row frow self.plague
        colors=np.empty(len(p),dtype="U10") #colors
        shapes=np.empty(len(q),dtype=matplotlib.path.Path) #shapes
        if len(p)!=len(q):
            print("ERROR IN COLORS_SHAPES")
            return -1
        for i in range(0,len(p)):
            if p[i]==1: #S
                colors[i]="tab:green"
            elif p[i]==-1: #E
                colors[i]="tab:orange"
            elif p[i]==-2: #C
                colors[i]="red"
            elif p[i]==-3: #D
                colors[i]="k"
            else: #R
                colors[i]="tab:blue"
                
            #does not look good but it works..
            if q[i]==0:
                shapes[i]=matplotlib.markers.MarkerStyle('o').get_path().transformed(
                    matplotlib.markers.MarkerStyle('o').get_transform())
            else:
                shapes[i]=matplotlib.markers.MarkerStyle('s').get_path().transformed(
                    matplotlib.markers.MarkerStyle('s').get_transform())
            
                
        return colors,shapes
        
        
    def plot_me(self,save=None,skip=1,interval=100,plot_seed=321,fsize=(16,9),msize=15,max_embedsize=100):
        #animated plotting
        frames=len(self.plague)
        fig = plt.figure(figsize=fsize)
        ax = fig.add_subplot(111)
        #plot with p0
        c,s=self.get_colors_shapes(self.plague[0],self.quarantine[0])
        sc=nx.draw_networkx_nodes(self.G,pos=nx.spring_layout(self.G,seed=plot_seed),node_size=msize,color=c)
        plt.axis("off")
        
        def draw(i):
            c,s=self.get_colors_shapes(self.plague[i],self.quarantine[i])
            sc.set_color(c)
            sc.set_paths(s)
            return sc,
        rc('animation', html='jshtml')
        rcParams['animation.embed_limit'] = max_embedsize #in MB
        anim=animation.FuncAnimation(fig, draw, frames=frames, interval=interval, blit=True)
        if save is not None:
            anim.save(save,writer="pillow")
        return anim
    
    def plot_init(self,plot_seed=321,fsize=(16,9),msize=15):
        #plot initial state without color code
        fig = plt.figure(figsize=fsize)
        ax = fig.add_subplot(111)
        nx.draw_networkx_nodes(self.G,pos=nx.spring_layout(self.G,seed=plot_seed),ax=ax,node_size=msize,color="tab:blue")
        
    def get_SIR_data(self):
        #get number of S,I,R ppl after each step
        S=np.array([])
        E=np.array([])
        R=np.array([])
        C=np.array([])
        D=np.array([])
        for j in range(0,len(self.plague)):
            s,i,r,c,d=self.sort_ppl(j)
            S=np.append(S,len(s))
            E=np.append(E,len(i))
            R=np.append(R,len(r))
            C=np.append(C,len(c))
            D=np.append(D,len(d))
        return S,E,R,C,D