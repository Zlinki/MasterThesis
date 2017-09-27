
#This notebook contains a quick implementation of the robust pca algorithm.
#We solve the problem using the implementation of "The Augmented Lagrange Multiplier
#method for exact recovery of corrupted low-rank matrices"
#In the end, we also use a solver to compare times and overall performance.
#We expect to find much more promising results using the implemantation on the latter paper:
#"/although interior point methods normally take very few iterations t converge, they have difficulty in
#handling large matrices because of the complexity of computing the step direction is O(m^6), where m is the dimension of 
#matrix [...] generic interior point solvers are too limited for Robust PCA ".

#Note that the robust PCA is used in many contexts. Ours is the clique problem which deals with nxn symetric matrices
#. However this implementation works for m x n matrices.

#Following the implementation of the solution proposed in the latter paper, we implement the solution
#proposed in "SHARP PERFORMANCE BOUNDS FOR GRAPH CLUSTERING VIA CONVEX OPTIMIZATION" by Korlakai, Oymak and Hassibi. for
#detection of cliques.




# The model
#D- a m x n matrix of data/variables. Usually D is nxn graph adjacency matrix.
# We seek to solve the problem min a||E||_1 + ||A||_* subject to E+A = D+I 
#a is a positive real number. 


#rpca is the main function. 
#If D is square matrix, we redefine D to be D = D+I.

function rpca(D,lambda,rho=1.6,stopCrit1=1.0e-7,stopCrit2=1.0e-5,maxIter=100)
    #Initialization
    #For problems where we have to add the identity...
    #dims =size(D)
    #if dims[1]==dims[2]
    #    D = D+I
    #end    
    mu = 1.25/norm(D)  
    E = zeros(D)
    k = 0
    Y = D/J(D,lambda)
    A=D
    #println("Comienzo!")
    #Iterations
    #In order to avoid evaluation of stoping criteria in the first loop, we do a while true and then 
    #add breaks if the stoping criteria are met.
    while true   
        # First solve A_k+1 = arg min L(A,E_k,Y_k,mu_k)
        singValDesc = svd(D-E+((1/mu)*Y))
        #SVD(A) returns a triple (U,S,V) where S contains the singular values and A=U*S*V' (' denotes ')
        #for SVD(A) to work correctly, A must be m x n where m>= n or else S wont have the proper dimensions.
        A = singValDesc[1]*perturb(diagm(singValDesc[2]),1/mu)*singValDesc[3]'
        #Now solve E_k+1= arg min L(A_k+1,E,Y_k,mu_k)
        perturbFactor = lambda*(1/mu)
        Eupdated = perturb(D-A+((1/mu)*Y),perturbFactor)
        Y = Y + mu*(D-A-Eupdated)
        mu=updateMu(mu,rho,E,Eupdated,D,stopCrit2)
        
       #Checks if both the first and second criterium are met.
        if  firstCriterium(stopCrit1,D,A,Eupdated) && secondCriterium(stopCrit2,D,mu,E,Eupdated) 
            println("Done")
            break
        end
        E=Eupdated
        k=k+1
        #println("Step $(k)")
        #println(mu)
        
        #forces the algorithm to stop if k>maxIter. This should never happen. maxIter is 1000 by default.
        if k>maxIter
            println("exceeded max number of iteration at solving")
            break
        end
    end
#    return(roundResultMatrix(A,mean(A)),roundResultMatrix(E,mean(E)))
     return(A,E)

end


#Stoping criteria functions  

function firstCriterium(stopCrit1,D,A,E)
    
    if vecnorm(D-A-E)/vecnorm(D)<stopCrit1
        return true
    end
    return false
end 

#Ek is the value of E computed on the kth step. Ek1 is the value of E computed in the k+1th step.
function secondCriterium(stopCrit2,D,mu,Ek,Ek1)
    
    if (mu*vecnorm(Ek1-Ek))/vecnorm(D)<stopCrit2
        return true
    end
    return false
end


    


# Updating functions
#Ek is the value of E computed on the kth step. Ek1 is the value of E computed in the k+1th step.
function updateMu(mu,rho,Ek,Ek1,D,epsilon2)
    
    eval = min(mu,sqrt(mu))*(vecnorm(Ek1-Ek) /vecnorm(D))
    if eval<epsilon2
        return rho*mu
    end
    return mu
end






#Other useful operators and methods

#perturbation operator. Perturbs every entry of the matrix W using the function f.
#W is a m x n matriz to perturb
#epsilon is the perturbation
function perturb(W,perturbation)
    map(W) do x
        f(x,perturbation)
    end
end


#J operator
#W is an m x n matriz
#lambda : parameter of the robust pca
function J(W,lambda)
    return max(norm(W,2),(1/lambda)*maximum(abs(W)))
end

#Method for perturbing x given epsilon. Used in the perturb function.
function f(x,epsilon)
    if x>epsilon 
        return x-epsilon
    end
    if x<-epsilon
        return x+epsilon
    else
        return(0)
    end
end

function roundResult(x,matrixMean)
    if abs(x)<matrixMean
        return 0
    end
    return 1
end

function roundResultMatrix(W,matrixMean)
    map(W) do x
        roundResult(x,matrixMean)
    end
end






#n is the dimension
#m is the sparcity
#k is the rank
function CreateRandomNoisy(n,m,k)
    X = randn(n,k)
    Y = randn(n,k)
    B = X*Y'
    totalEntries = n*n
    randomEntries1=rand(1:n,n*n-m)    
    randomEntries2=rand(1:n,n*n-m)
    A= randn(n,n)
    randomEntries1=rand(1:n,m)
    randomEntries2=rand(1:n,m)
    S = zeros(A)
    for i = 1:m
        S[randomEntries1[i],randomEntries2[i]]=A[randomEntries1[i],randomEntries2[i]]   
    end
    return (B,S)
end 

#Function that test if the recover was succsesfull given the known matrices.
function testMethod(n,m,k,epsilon)
    for i in 0.0:0.1:0.9
    println(i)
        j=i/(1-i)
        real= CreateRandomNoisy(n,m,k)
        solve=rpca(real[1]+real[2],j,1.6,1.0e-7,1.0e-5,5000)
        if checkSolution(real[1],real[2],solve[1],solve[2],epsilon)==true
            println("SolucionEncontrada!")
            return(real[1],real[2],solve[1],solve[2])
        end
    end
        println("No encontre solucion :/")
end

#Function that checks if a given pair of matrices is a solution
function checkSolution(RealLowRank,RealSparse,LowRank,Sparse,epsilon)
    tolerance =  vecnorm(LowRank-RealLowRank)/vecnorm(RealLowRank)+vecnorm(Sparse-RealSparse)/vecnorm(RealSparse)
    if tolerance < epsilon
        return true
    end
        return false
end    

    
    
    

testMethod(10,3,2,0.1)


