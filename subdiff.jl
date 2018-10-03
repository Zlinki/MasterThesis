A=[1 1 -1 -1 -1 ; 1 1 -1 -1 -1 ; -1 -1 1 1 1 ; -1 -1 1 1 1 ; -1 -1 1 1 1 ]
display(A)
E= eigvals(A)
D, V = eig(A)
D2 = Diagonal(D)
display(V*D2*transpose(V))
R=Diagonal([-1,-1,-1,-1,1])
display(R)
S=V*R*transpose(V)
display(S)
display(round(S))
trace(S*A)
print(eigvals(S))
display(V)
display(S)#Queremos S negativo en el cluster y positivo afuera

#More examples should be good.
u = [1 1 1 1 1 1 -1 -1 -1 -1 -1]
A=transpose(u)*u
display(A)
D, V = eig(A)
D2 = Diagonal(D)
display(V*D2*transpose(V))
R=Diagonal([-1,1,1,1,1,1,1,1,1,1,1])
S=V*R*transpose(V)
display(round(S,1))#This one seems to have the correct signs!
display(S)#Conjecture: The matrix R with eigenvalues given by
trace(A*S)
norm(S)

A=[ 1 2 0 ; 1 1 1 ; 3 -4 2]
b=transpose([1 1 3])
display(A^(-1)*b)

#Looking for orthogonal exchanges:
n=7

function S(i,j)
  u =[]
  for t in 1:n
    if t==i append!(u,1)
    elseif t==j append!(u,-1)
    else
      append!(u,0)
    end
  end
  return u*transpose(u)
end




display(S(1,2))

v= [1 1 1 1 1 1 1]
G = transpose(v)*v
N=S(1,2)-S(1,3)
display(N)
N2 = S(1,4)-S(1,5)
N3 = S(1,6)-S(1,7)
display(G*N)
display(N*G)

trace(N3)

norm(N-N2-N3)


display(N3*N)

norm(N+N2)

max(norm(N),norm(N2))

norm(S(1,2)-S(1,3))

norm((S(1,2)-(S(1,3)+S(1,4)+S(1,5)+ S(1,6))/4))
norm(S(1,2)-(S(1,3)+S(1,4))/2)

norm(S(1,2)-(S(1,3)))
