export write_full_data, to_crys, unfold_mmn


function to_crys(nnkp, bvectors)
    bvectors_cryst = [zeros(3) for _ in 1:size(bvectors)[1]]
    for ib in 1:size(bvectors)[1]
        bvectors_cryst[ib] = inv(nnkp.recip_lattice)*bvectors[ib]
    end
    return bvectors_cryst
end 



function write_full_data(seedname)
    nnkp = read_nnkp(seedname*".nnkp")
    Kstencil = read_nnkp_compute_bweights(seedname*".nnkp")
    isym = read_w90_isym(seedname*".isym")
    rescale_repmat(isym.gk) #not sure why this is needed (to check)
    f2i = fbz2ibz(nnkp.kpoints, isym.kpoints, isym.symops)
    println(f2i)
    #Eig
    #eig = read_eig(seedname*".ieig")
    #write_eig(f2i, eig, seedname*".eig")

    #Amn
    A  = read_amn(seedname*".iamn")
    #fill vector for the factor e^{-ikR_n'} appearing in Koretsune (9) 
    Rn = find_pos(isym.symops,  
                nnkp.projections,
                isym.D
    )
    println(Rn)
    Asymm = symmetrize(A, isym, Rn)
    Afull = unfold_amn(Asymm, f2i, Rn, isym)

    #Mmn
    M, kpb_k, kpb_G = read_mmn(seedname*".immn")
    bvectors_cryst = to_crys(nnkp, Kstencil.bvectors)
    Mfull = unfold_mmn(M, kpb_k, kpb_G, nnkp, isym, f2i)


####
#####check Amn with symwannier
####Apy = read_amn(seedname*".amn")
####for ik in 1:size(f2i)[1]
####    i = findall( >(1.e-07), abs.(Afull[ik]-Apy[ik]))
####    println("ik = ", ik , "    " , i)
####end 

end




function unfold_mmn(M, kpb_k, kpb_G, nnkp, isym, f2i)
    Mfull = [[zeros(ComplexF64, size(M[1][1])) for _ in 1:size(M[1])[1]] for _ in 1:size(M)[1]]

    #for ikf in 1:size(nnkp.kpoints)[1]
    #    kf       = nnkp.kpoints[ikf]
    #    iki      = f2i[ikf][1][1]
    #    isym_ikf = f2i[ikf][1][2]
    #    ki       = isym.kpoints[iki]
    #    
    #    for ibf in 1:size(bvectors_cryst)[1]
    #        # getting b point as (k+b) + G - k 
    #        ikbf = nnkp.kpb_k[ikf][ibf] #index of k+b
    #        G = nnkp.kpb_G[ikf][ibf]
    #        bf  = nnkp.kpoints[ikpb] + G - nnkp.kpoints[ikf] # unrotated b vector
    #        #bi = h^{-1} bf
    #        bi = rotate(bf, isym.symops[isym.symops[isym_ikf].invs]) 
    #        ibi = find(bi, bvectors_cryst)
    #        ikbi= find( ki + bi, nnkp.kpoints)
    #        
    #        ikbi_ibz  = f2i[ikbi][1][1]            
    #        isym_ikbi = f2i[ikbi][1][2] #symmetry to get ki+bi in IBZ
    #        isym_ikbf = f2i[ikbf][1][2] #symmetry to get kf+bf in IBZ
    #        #get equivalent operation of g^{-1}[isym_ikbi]*g^{-1}[isym_ikf]g[isym_ikbf] 
    #        isym_h, tDelta = get_symop(isym.symops, [[isym_ikbi, -1], [isym_ikf, -1], [isym_ikbf, +1]])
    #        println(ikf, " ", bf, " ", bi, " ", isym_ikf, " " ,isym.symops[isym_ikf].t_rev," ",isym.symops[isym_ikf].invs )#isym_tot, " ", tDelta)
#
    #        factor = 1#changes with spinors
    #        if isym.symops[isym_ikbi].t_rev == 1
    #            Mfull[ikf][ibf] = M[iki][ibi]*conj.(isym.g[ikbi_ibz].h[isym_h])*factor
    #        else
    #            Mfull[ikf][ibf] = M[iki][ibi]*isym.g[ikbi_ibz].h[isym_h]*factor
    #        end
#
    #        if isym.symops[isym_ikf].t_rev == 1
    #            Mfull[ikf][ibf] = conj.[Mfull[ikf][ibf]]
    #        end
#
    #        kbi_ibz = isym.kpoints[ikbi_ibz]
    #        Mmn[ikf][ibf] = Mmn[ikf][ibf]*exp(-im*2.0*π*(dot(bi, isym.symops[isym_ikf].ft) + dot(kbi_ibz, tDelta)))
    #    end
    #end
end