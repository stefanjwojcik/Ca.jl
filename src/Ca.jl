module Ca

using Pkg, CSV, LinearAlgebra, DataFrames, DataFramesMeta, InvertedIndices

# methods to export and make available
export foo, bar, ca

foo(x::T, y::T) where T <: Real = x + y - 5
bar(z::Float64) = foo(sqrt(z), z)

# the initial correspondence analysis function
function ca(M::Array{Float64, 2}, k::Integer)
    O = M ./ sum(M)
    #row totals (row masses) by each of the column totals.
    Ep = sum(O, dims=2) .* sum(O, dims=1)
    #En = Ep * sum(M)
    Z = (O .- Ep) ./ sqrt.(Ep)
    U, S, V = svd(Z)
    U[:, k]
end

function ca(M::Array{Int64, 2}, k::Integer)
    # convert if not a float matrix
    M = convert(Array{Float64, 2}, M)
    O = M ./ sum(M)
    #row totals (row masses) by each of the column totals.
    E = sum(O, dims=2) .* sum(O, dims=1)
    Z = (O .- E) ./ sqrt.(E)
    U, S, V = svd(Z)
    U[:, k]
end

function ca(M, k)
    println("Boo. Your args are not matrices or integers, I cannot compute 💻")
end


include("data/cadata.jl")


###### REVISING THE CA FUNCTION

function Ca_V2 (obj::Array{Int64, 2};
    nd::Integer = missing,
    cnames = missing,
    rnames = missing,
    suprow = missing,
    supcol = missing,
    subsetrow = missing,
    subsetcol = missing)

    # assign things
    nd0 = nd
    I = size(obj)[1]
    J = size(obj)[2]
    rn = rnames
    cn = cnames
    #N <- matrix(as.matrix(obj), nrow = I, ncol = J) already a matrix
    N = obj
    Ntemp = N
    # Temporary holders based on supplementary rows/cols
    NtempC = NtempR = N
    suprow = sort(suprow)
    supcol = sort(supcol)
    if !ismissing(supcol[1]) & !ismissing(suprow[1])
         # NEED TO FIX NOT PHRASING
         NtempC = NtempC[Not(suprow), :)
         NtempR = Ntemp[: , Not(supcol)]
    end
    if !ismissing(supcol[1])
         SC = Matrix(NtempC[:, supcol])
         Ntemp = Ntemp[:, Not(supcol)]
         cs.sum = apply(SC, 2, sum)
    end

    if !ismissing(suprow[1])
         #SR <- matrix(as.matrix(NtempR[suprow, ]), nrow = length(suprow))
         Ntemp = Ntemp[-suprow, ]
         rs.sum = apply(SR, 1, sum)
    end
    #N <- matrix(as.matrix(Ntemp), nrow = dim(Ntemp)[1], ncol = dim(Ntemp)[2])
    subsetrowt = subsetrow
    if !ismissing(subsetrow[1]) & !ismissing(suprow[1])
         subsetrowi = subsetrow
         subsetrowt = sort(vcat(subsetrow, suprow))
         subsetrowt = subsetrowt[Int.(Set(subsetrowt))]
         I = length(subsetrowt)
         for q in reverse(1:length(suprow))
             subsetrow = subsetrow[subsetrow ∉ suprow[q]]
             subsetrow = subsetrow - 1*(suprow[q] < subsetrow)
         end
        suprow = [q for q in 1:length(suprow) if suprow[q] ∈ subsetrowt]
     end

     subsetcolt = subsetcol

     if !ismissing(subsetcol[1]) & !ismissing(supcol[1])
         subsetcoli = subsetcol
         subsetcolt = sort(vcat(subsetcol, supcol))
         subsetcolt = subsetcolt[Int.(Set(subsetcolt))]
         J = length(subsetcolt)
     for q in reverse(1:length(supcol))
         subsetcol = subsetcol[subsetcol ∉ supcol[q]]
         subsetcol = subsetcol - 1*(supcol[q] < subsetcol)
     }
     for q in 1:length(supcol)
         supcol = [q for q in 1:length(subsetcol) if supcol[q] ∈ subsetcolt]
         ## LEFT OFF HERE
     dim.N = size(N)
     if !ismissing(subsetrow[1])) {
     if (!is.na(supcol[1]))
       SC <- as.matrix(SC[subsetrow, ])
    }
    if (!is.na(subsetcol[1])) {
     if (!is.na(suprow[1]))
       SR <- matrix(as.matrix(SR[, subsetcol]), nrow = length(suprow))
    }
    if (is.na(subsetrow[1]) & is.na(subsetcol[1])) {
    nd.max <- min(dim.N) - 1
    }
    else {
    N00 <- N
    if (!is.na(subsetrow[1]))
      N00 <- N00[subsetrow, ]
    if (!is.na(subsetcol[1]))
      N00 <- N00[, subsetcol]
    dim.N <- dim(N00)
    nd.max <- min(dim.N)
    if (!is.na(subsetrow[1]) & is.na(subsetcol[1])) {
      if (dim.N[1] > dim.N[2])
        nd.max <- min(dim.N) - 1
    }
    else {
      if (is.na(subsetrow[1]) & !is.na(subsetcol[1])) {
        if (dim.N[2] > dim.N[1]) {
          nd.max <- min(dim.N) - 1
        }
      }
    }
    }
  if (is.na(nd) | nd > nd.max)
    nd <- nd.max
  n <- sum(N)
  P <- N/n
  rm <- apply(P, 1, sum)   row masses
  cm <- apply(P, 2, sum)   column masses
  eP <- rm %*% t(cm)   expected probability - rc times cm
  eN <- eP * n
  S <- (P - eP)/sqrt(eP)
  if (!is.na(subsetcol[1])) {
    S <- S[, subsetcol]
    cm <- cm[subsetcol]
    cn <- cn[subsetcolt]
  }
  if (!is.na(subsetrow[1])) {
    S <- S[subsetrow, ]
    rm <- rm[subsetrow]
    rn <- rn[subsetrowt]
  }
  chimat <- S^2 * n
  dec <- svd(S)
  sv <- dec$d[1:nd.max] # eigenvalues of SVM
  ev <- sv^2
  cumev <- cumsum(ev)
  totin <- sum(ev)

  u <- dec$u
  v <- dec$v

  rin <- apply(S^2, 1, sum)
  cin <- apply(S^2, 2, sum)
  rachidist <- sqrt(rin/rm)
  cachidist <- sqrt(cin/cm)
  rchidist <- rep(NA, I)
  cchidist <- rep(NA, J)
  if (!is.na(subsetrow[1])) {
    obj <- obj[subsetrowt, ]
  }
  if (!is.na(subsetcol[1])) {
    obj <- obj[, subsetcolt]
  }
  if (!is.na(suprow[1])) {
    if (is.na(supcol[1])) {
      P.stemp <- matrix(as.matrix(obj[suprow, ]), nrow = length(suprow))
    }
    else {
      P.stemp <- matrix(as.matrix(obj[suprow, -supcol]),
                        nrow = length(suprow))
    }
    P.stemp <- P.stemp/apply(P.stemp, 1, sum)
    P.stemp <- t((t(P.stemp) - cm)/sqrt(cm))
    rschidist <- sqrt(apply(P.stemp^2, 1, sum))
    rchidist[-suprow] <- rachidist
    rchidist[suprow] <- rschidist
  }
  else rchidist <- rachidist
  if (!is.na(supcol[1])) {
    if (is.na(suprow[1])) {
      P.stemp <- as.matrix(obj[, supcol])
    }
    else {
      P.stemp <- as.matrix(obj[-suprow, supcol])
    }
    P.stemp <- t(t(P.stemp)/apply(P.stemp, 2, sum))
    P.stemp <- (P.stemp - rm)/sqrt(rm)
    cschidist <- sqrt(apply(P.stemp^2, 2, sum))
    cchidist[-supcol] <- cachidist
    cchidist[supcol] <- cschidist
  }
  else {
    cchidist <- cachidist
  }
  phi <- as.matrix(u[, 1:nd])/sqrt(rm)
  gam <- as.matrix(v[, 1:nd])/sqrt(cm)
  if (!is.na(suprow[1])) {
    cs <- cm
    gam.00 <- gam
    base2 <- SR/matrix(rs.sum, nrow = nrow(SR), ncol = ncol(SR))
    base2 <- t(base2)
    cs.0 <- matrix(cs, nrow = nrow(base2), ncol = ncol(base2))
    svphi <- matrix(sv[1:nd], nrow = length(suprow), ncol = nd,
                    byrow = TRUE)
    base2 <- base2 - cs.0
    phi2 <- (t(as.matrix(base2)) %*% gam.00)/svphi
    phi3 <- matrix(NA, ncol = nd, nrow = I)
    phi3[suprow, ] <- phi2
    phi3[-suprow, ] <- phi
    rm0 <- rep(NA, I)
    rm0[-suprow] <- rm
    P.star <- SR/n
    rm0[suprow] <- NA
    rin0 <- rep(NA, I)
    rin0[-suprow] <- rin
    rin <- rin0
    rm.old <- rm
    rm <- rm0
  }
  if (!is.na(supcol[1])) {
    if (!is.na(suprow[1])) {
      rs <- rm.old
    }
    else {
      rs <- rm
    }
    phi.00 <- phi
    base2 <- SC/matrix(cs.sum, nrow = nrow(SC), ncol = ncol(SC),
                       byrow = TRUE)
    rs.0 <- matrix(rs, nrow = nrow(base2), ncol = ncol(base2))
    svgam <- matrix(sv[1:nd], nrow = length(supcol), ncol = nd,
                    byrow = TRUE)
    base2 <- base2 - rs.0
    gam2 <- (as.matrix(t(base2)) %*% phi.00)/svgam
    gam3 <- matrix(NA, ncol = nd, nrow = J)
    gam3[supcol, ] <- gam2
    gam3[-supcol, ] <- gam
    cm0 <- rep(NA, J)
    cm0[-supcol] <- cm
    P.star <- SC/n
    cm0[supcol] <- NA
    cin0 <- rep(NA, J)
    cin0[-supcol] <- cin
    cin <- cin0
    cm <- cm0
  }
  if (exists("phi3")) {
    phi <- phi3
  }
  if (exists("gam3")) {
    gam <- gam3
  }
  dims <- paste0("Dim", seq_along(sv))[1:nd]
  dimnames(phi) <- list(rn, dims)
  dimnames(gam) <- list(cn, dims)
  ca.output <- list(sv = sv, nd = nd0, rownames = rn, rowmass = rm,
                    rowdist = rchidist, rowinertia = rin, rowcoord = phi,
                    rowsup = suprow, colnames = cn, colmass = cm, coldist = cchidist,
                    colinertia = cin, colcoord = gam, colsup = supcol, N = N,
                    call = match.call())
  class(ca.output) <- "ca"
  return(ca.output)
}

end # module
