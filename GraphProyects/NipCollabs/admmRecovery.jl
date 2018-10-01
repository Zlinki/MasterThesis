
#This notebook contains an implementation to solve the problem min Sum ||A-B_i||_1 + ||2A-11^t||_nucl.
#We solve the problem an the implementation based on "The Augmented Lagrange Multiplier
#method for exact recovery of corrupted low-rank matrices" and some notes of Boyd for ADMM.
#which we link here. https://web.stanford.edu/~boyd/papers/pdf/admm_slides.pdf

#"/although interior point methods normally take very few iterations t converge, they have difficulty in
#handling large matrices because of the complexity of computing the step direction is O(m^6), where m is the dimension of 
#matrix [...] generic interior point solvers are too limited for Robust PCA ".

#Note that the robust PCA is used in many contexts. Ours is the clique problem which deals with nxn symetric matrices
#. However this implementation works for m x n matrices. Our implementation is an extension of the robust pca algorithm.


#admmRecovery does not call mosek. It uses the prox operation to find xi+1 and Z.
function admmRecovery(observations,lambda,n,rho=1.6,stopCrit1=1.0e-5,maxIter=200)
    
    #Important! all matrices in the observations are assumed to have 1 in the diagonal.
    #Transformation of the B_i to the C_i where C_i = 2B_i-11^t
    observationsTransformed = transform(observations,n)
    #Initialization (All parameters are based on the paper.)
    k = 0
    N = length(observationsTransformed)
    averObservations=AvergeMatrix(observationsTransformed)
    Y = Array{Any}(N)
    [Y[r]=observationsTransformed[r]/J(observationsTransformed[r],lambda) for r in 1:N]
    Z= averObservations
    softZ = zeros(Z)
    Zviejo= zeros(Z)
    Xviejo = Array{Any}(N)
    [Xviejo[r]=zeros(Z) for r in 1:N]
    
    Xnuevo = Array{Any}(N)
    [Xnuevo[r]=zeros(Z) for r in 1:N]
    
    #mu = 1.25/norm(averObservations,2)
    #Iterations
    #In order to avoid evaluation of stoping criteria in the first loop, we do a while true and then 
    #add breaks if the stoping criteria are met.
    while true 
                
        # First solve x^i_k+1 = arg min 1/2||x^i_k-C_i||+mu/2||x^i_k-Z_k+y^i_k||_2^2
        println("solving for x...")
        [Xnuevo[i] = solveForXprox(observationsTransformed[i],Y[i],Z,rho,n) for i in 1:N]
       

        println("Solved, now solving for Z")
        tic()
         # now we solve Z_k+1 = arg min prox_g,rho (1/N sum(x_i_k+1 + 1/rho y_i^k))
        singValDesc = svdfact!(AvergeMatrix(Xnuevo +(1/rho)*Y))
        #SVD(A) returns a triple (U,S,V) where S contains the singular values and A=U*S*V' (' denotes ')
        #for SVD(A) to work correctly, A must be m x n where m>= n or else S wont have the proper dimensions.
        softZ = diagm(SoftThreshold.(singValDesc[:S],lambda/rho))
        softZ = convert(Array{Float64,2},softZ)    
        Z = singValDesc[:U]*softZ*singValDesc[:Vt]       
        toc()
        println("solved. now solving Y")
        #update Y
        [Y[i] = Y[i]+ rho*(Xnuevo[i]-Z)for i in 1:N]
        
        
        #compute the residuals to check stopping conditions and compute rho
        residual = AvergeMatrix(Xnuevo)-Z
        residualDual = rho*(Zviejo-Z)
        #check stopping condition
     #   println("Cheking stop criteriums and updating.")
     #   if  firstCriterium(stopCrit1,n,residual,residualDual)
     #       println("Done")
     #       break
     #   end
        #update rho
        rho = updaterho(rho,residual,residualDual)
         
        #update Xviejo
        Xviejo = Xnuevo
        #update Zviejo
        Zviejo = Z
        # update k 
        println("Step $(k)")
    #    println(norm(Xnuevo[1],2))
        k = k+1
        
        println(rho)
        if k>maxIter
            println("exceeded max number of iteration at solving")
            break
        end    
    end
    return((Z,Xnuevo))
end
        
        
        

function solveForXprox(obs,Y,Z,mu,n)
    actualizado = SoftThreshold(Z-(1/mu)*Y-obs,1/(2*mu))+obs
    return(actualizado)
end    


#Aux functions


#perturbation operator. Perturbs every entry of the matrix W using the function f.
#W is a m x n matriz to perturb
#epsilon is the perturbation
function SoftThreshold(W,perturbation)
    map(W) do x
        perturb(x,perturbation)
    end
end



#Method for perturbing x given epsilon. Used in the softThresholding function.
function perturb(x,epsilon)
    if x>epsilon 
        return x-epsilon
    end
    if x<-epsilon
        return x+epsilon
    else
        return(0)
    end
end



#Method to average matrices

function AvergeMatrix(observations)
    size = length(observations)
    return (1/size)*sum(observations)
end

#J operator
#W is an m x n matriz
#lambda : parameter of the algorithm
function J(W,lambda)
    return max(norm(W,2),(1/lambda)*maximum(abs.(W)))
end

#Function  transform the observed matrices to the desired form. All observed matrices
# have 1 on the diagonal

function transformMatrix(mat,n)
    return(2*mat-ones(n)*ones(n)') 
end

function transform(observations,n)
    return(transformMatrix.(observations,n))
end



function firstCriterium(stopCrit1,n,residual,residualDual)
    
    if  sqrt(n)*norm(residual,2)<= stopCrit1 && sqrt(n)*norm(residualDual,2)<= stopCrit1
        return(true)
    end
    return false
end

# Updating functions

function updaterho(rho,residual,residualDual)
    rhoNew = rho
    if norm(residual,2)>=10*norm(residualDual,2)
    rhoNew = 2*rho
    end
    if norm(residualDual,2)>=10*norm(residual,2)
        rhoNew = rho/2
    end
    return(rhoNew)
end
        


#using NBInclude
#nbinclude("AuxFunctions.ipynb")
#nbinclude("GraphLearning.ipynb")

#srand(12536)


#grafo20obs5 = Array{Any}(5)
#grafo20obs5[1]=erdosgraph([0.77,0.77],0.1,[12,8],20)+diagm(ones(20))
#grafo20obs5[2]=erdosgraph([0.77,0.77],0.1,[12,8],20)+diagm(ones(20))
#grafo20obs5[3]=erdosgraph([0.77,0.77],0.1,[12,8],20)+diagm(ones(20))
#grafo20obs5[4]=erdosgraph([0.77,0.77],0.1,[12,8],20)+diagm(ones(20))
#grafo20obs5[5]=erdosgraph([0.77,0.77],0.1,[12,8],20)+diagm(ones(20))
#1

#0.28889 s
#tic()
#recuperaRegular = learnRegularized2A([1/5,1/5,1/5,1/5,1/5],1.6,grafo20obs5,20)
#toc()

#recuperaRegular[1]

#function admmRecovery2(observations,lambda,n,rho=1.6,stopCrit1=1.0e-7,stopCrit2=1.0e-5,maxIter=10000)


#64.3 segundos
#tic()
#respuestaAdmm = admmRecovery2(grafo20obs5,1.6,20,1.6,1.0e-5,1.0e-5,500)
#toc()

#bli =(respuestaAdmm[1]+ones(20)*ones(20)')/2

#0.069 s
#tic()
#respuestaAdmmProx = admmRecovery(grafo20obs5,1.6,20,1.6,1.0e-5,500)
#toc()

#bli =(respuestaAdmmProx[1]+ones(20)*ones(20)')/2

#grafo200obs3 = Array{Any}(3)
#grafo200obs3[1]=erdosgraph([0.85,0.85],0.1,[100,100],200)+diagm(ones(200))
#grafo200obs3[2]=erdosgraph([0.85,0.85],0.1,[100,100],200)+diagm(ones(200))
#grafo200obs3[3]=erdosgraph([0.85,0.85],0.1,[100,100],200)+diagm(ones(200))


#7 horas
#tic()
#respuestaAdmmGrande = admmRecovery2(grafo200obs3,1.6,200,1.6,1.0e-5,1.0e-5,200)
#toc()

#(respuestaAdmmGrande[1]+ones(200)*ones(200)')/2

#3.15 s (antes de la correccion del cast tardabamos 7 minutos)
# sin embargo notamos que desde el step 20 ya no cambia la norma de la matriz asi que volvemos a tratar mas abajo...
#tic()
#respuestaAdmmProx = admmRecovery(grafo200obs3,1.6,200,1.6,1.0e-5,200)
#toc()

#(respuestaAdmmProx[1]+ones(200)*ones(200)')/2

#grafo1000obs3 = Array{Any}(3)
#grafo1000obs3[1]=erdosgraph([0.85,0.85],0.1,[250,750],1000)+diagm(ones(1000))
#grafo1000obs3[2]=erdosgraph([0.85,0.85],0.1,[250,750],1000)+diagm(ones(1000))
#grafo1000obs3[3]=erdosgraph([0.85,0.85],0.1,[250,750],1000)+diagm(ones(1000))


#139 min = 2h
#tic()
#respuestaAdmmProx = admmRecovery(grafo1000obs3,1.6,1000,1.6,1.0e-5,1.0e-5,30)
#toc()

#bli =(respuestaAdmmProx[1]+ones(1000)*ones(1000)')/2

#corrigiendo con el cast
#9.26 min
#tic()
#respuestaAdmmProx = admmRecovery(grafo1000obs3,2,1000,1.6,1.0e-5,200)
#toc()


#bli =(respuestaAdmmProx[1]+ones(1000)*ones(1000)')/2


