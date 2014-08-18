using Devectorize

function async_gs(m,b,x,r,eps)
###
# Solves the system m * x = b using distributed async Gauss-Seidel method
# m is a distributed matrix
# b,x,r are distributed vectors
# eps is the stopping criterium
###
    res = Inf
    I = localindexes(x)[1]
    localx = convert(Array,x)
    n = size(m,2)

    updatingx = false
    calculatingr = false

    iter = 0

    while res > eps

        # Updating x

        for i in I
            s = 0
            for j = 1:n
                if i != j
                    s += m[i,j] * localx[j]
                end
            end

            localx[i] = (b[i] - s) / m[i,i]
            localpart(x)[i - I.start + 1] = localx[i]
        end

        # Stopping criteria

        res = 0

        for i in I

            s = b[i]

            for j = 1:n

                s -= m[i,j] * localx[j]

            end

            localpart(r)[i - I.start + 1] = abs(s)

            res = max(res,s)

        end

        if res <= eps && !calculatingr

            calculatingr = true
            @async begin
                res = norm(r,Inf)
                calculatingr = false
            end                

        end

        if !updatingx

            updatingx = true
            @async begin
                copy!(localx,x)
                updatingx = false
            end

        end

        iter += 1
    end
    return iter
end

function shared_gs(m,b,x,r,eps)
    ###
    # Solves the system m * x = b using a shared memory parallel
    # Gauss-Seidel method
    ###
    res = norm(r,Inf)
    I = localindexes(x)
    n = size(m,2)

    iter = 0
    while res > eps

        # Updating x

        for i in I
            s = 0
            for j = 1:n
                if i != j
                    s += m[i,j] * x[j]
                end
            end

            x[i] = (b[i] - s) / m[i,i]

        end

        # Stopping criterium

        res = 0

        for i in I

            s = b[i]

            for j = 1:n

                s -= m[i,j] * x[j]

            end

            r[i] = abs(s)

            res = max(res,s)

        end

        if res <= eps

            res = norm(r,Inf)

        end

        sleep(rand() / 1000)

        iter += 1
    end
    return iter
end

function serial_gs(m,b,x,r,eps)

    res = Inf
    n = size(m,1)
    iter = 0

    while res > eps
        for i in 1:n
            s = b[i]
            for j = 1:n
                if i != j
                    s -= m[i,j] * x[j]
                end
            end
            x[i] = s / m[i,i]
        end

        res = 0

        for i = 1:n
            s = b[i] - BLAS.dot(n,m[i,:],1,x,1)
#            for j = 1:n
#                s -= m[i,j] * x[j]
#            end
            res = max(res,abs(s))
        end

        iter += 1
    end
    
    return iter
end

function gs(m,b,x,N)

    for k = 1:N
        I = localindexes(x)[1]
        for i in I
            s = 0
            for j = 1:size(m,2)
                if i != j
                    s += m[i,j] * x[j]
                end
            end
            localpart(x)[i - I.start + 1] = (b[i] - s) / m[i,i]
        end
    end
end

function initializeRandomDistributedDDMatrix(I,l,u)
    m = l + (u - l) * rand(size(I[1],1),size(I[2],1))
    for i in 1:size(I[1],1)
        @devec m[i,I[1][i]] = sum(abs(m[i,:]))
    end
    return m
end

function generateRandomDistributedSystem(n,sol=ones((n,1)),w=workers(),
                                         l=-100,u=100)
    
    dist = [size(w,1),1]

    M = DArray((I)->initializeRandomDistributedDDMatrix(I,l,u),(n,n),w,dist);

    b = dzeros((size(M,1),1),w,dist)
    r = similar(b)

    for i in w
        @spawnat i localpart(b)[:,1] = localpart(M) * sol
        @spawnat i localpart(r)[:,1] = b[localindexes(b)[1],1]
    end

    x0 = dzeros((size(M,1),1),w,dist)

    return (M,b,x0,r)
end

# Para rodar
# require("gs_distribuido.jl")
#
# (M,b,x,r) = generateRandomDistributedSystem(20)
#
# for i in workers()
#     @spawnat i async_gs(M,b,x,r,0.1)
# end
#
# Ou, a versao de 1 linha (imprime o tempo e o total de iteracoes):
#
# @time reduce(+,pmap(fetch,{(@spawnat i async_gs(M,b,x,r,0.1)) for i in workers()}))
