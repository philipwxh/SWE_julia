function swe_2d_rhs(U,ops,dis_cst,vgeo,fgeo,nodemaps, dt, tol, g, pre_allo)
    # unpack args
    h, hu, hv, btm = U
    (Qr_ID,Qs_ID,Qrb_ID,Qsb_ID,Qr_ES,Qs_ES,QNr_sbp, QNs_sbp,
     E,M_inv, Mf_inv, Pf,Fmask ) = ops
    Cf, C, C_x, C_y, Nq, Nfq, nelem= dis_cst
    rxJ,sxJ,ryJ,syJ,J = vgeo
    nxJ,nyJ,sJ,nx,ny = fgeo
    mapP,mapB = nodemaps
    (u, v, uf, vf, hf, huf, hvf, hP, huP, hvP, dh, dhu, dhv,
    rhs1_ID, Mrhs1_ID, rhs1_ES, Mrhs1_ES, lf, lff,
    f1_ES, f2_ES, f3_ES, f1_ID, f2_ID, f3_ID, rhs1_CL, rhs2_CL, rhs3_CL,
    lambdaf, lambdaP, cf, UL_E, UR_E, vgeo_e,
    h_L_next, h_L_next_f, f1_IDf, f1_ESf, FS1, FS2, FS3)= pre_allo
    co_d = 4
    if dt == 0
        fill!(rhs1_CL, zero(eltype(rhs1_CL)));
        fill!(rhs2_CL, zero(eltype(rhs2_CL)));
        fill!(rhs3_CL, zero(eltype(rhs3_CL)));
        return rhs1_CL, rhs2_CL, rhs3_CL;
    end
    # @batch parallelizes
    @batch for e = 1:nelem
        for i = 1:Nq
            u[i,e] = hu[i,e]/h[i,e];
            v[i,e] = hv[i,e]/h[i,e];
        end
        for i = 1:Nfq
            hf[i,e]  = h[Fmask[i], e];
            huf[i,e] = hu[Fmask[i], e];
            hvf[i,e] = hv[Fmask[i], e];
            # uf[i,e] = huf[i,e]/hf[i,e];
            # vf[i,e] = hvf[i,e]/hf[i,e];
            lambdaf[i,e] = abs(huf[i,e]/hf[i,e]*nx[i,e]+
                   hvf[i,e]/hf[i,e]*ny[i,e]) + sqrt(g*hf[i,e])
        end
    end
    @batch for e = 1:nelem
        for i = 1:Nfq
            hP[i,e] = hf[mapP[i,e]];
            huP[i,e] = huf[mapP[i,e]];
            hvP[i,e] = hvf[mapP[i,e]];
        end
    end
    @batch for i in mapB
        vu_dot_n  = huf[i]*nx[i]+ hvf[i]*ny[i];
        vu_dot_n = vu_dot_n;
        huP[i] =  huf[i] - 2*vu_dot_n*nx[i];
        hvP[i] =  hvf[i] - 2*vu_dot_n*ny[i];
        # huP[i] =  -huf[i];
        # hvP[i] =  -hvf[i];
    end
    @batch for e = 1:nelem
        for i = 1:Nfq
            dh[i,e] = hP[i,e]-hf[i,e];
            dhu[i,e] = huP[i,e]-huf[i,e];
            dhv[i,e] = hvP[i,e]-hvf[i,e];
            lambdaP[i,e] = lambdaf[mapP[i,e]];
            cf[i,e] = max(lambdaP[i,e], lambdaf[i,e]);
        end
    end
    ##surface part
    #ES part
    @batch for e = 1:nelem
        for i = 1:Nfq
            ul = (hf[i,e], huf[i,e], hvf[i,e]);
            ur = (hP[i,e], huP[i,e], hvP[i,e]);
            du = (dh[i,e], dhu[i,e], dhv[i,e]);
            f_geo = (nxJ[i,e], nyJ[i,e], sJ[i,e], cf[i,e]);
            FS1[i,e], FS2[i,e], FS3[i,e] = swe_2d_esdg_n_surface(ul, ur, du, f_geo, g);
        end
    end

    # fill!(f1_ES, zero(eltype(f1_ES)));
    # fill!(f2_ES, zero(eltype(f2_ES)));
    # fill!(f3_ES, zero(eltype(f3_ES)));
    # @batch for e = 1:nelem
    #     for i = 1:Nfq # loop over surface nodes
    #         j = Fmask[i] # face node corresponding to volume node
    #         f1_ES[i, e] += wf[j] * FS1[j, e]
    #         f2_ES[i, e] += wf[j] * FS2[j, e]
    #         f3_ES[i, e] += wf[j] * FS3[j, e]
    #     end
    # end
    mul!(f1_ES, Pf, FS1);mul!(f2_ES, Pf, FS2); mul!(f3_ES, Pf, FS3);

    # fill!(f1_ID, zero(eltype(f1_ES)));
    # fill!(f2_ID, zero(eltype(f2_ES)));
    # fill!(f3_ID, zero(eltype(f3_ES)));
    #ID part
    @batch for e = 1:nelem
        for i = 1:Nfq
            ul = (hf[i,e], huf[i,e], hvf[i,e]);
            ur = (hP[i,e], huP[i,e], hvP[i,e]);
            du = (dh[i,e], dhu[i,e], dhv[i,e]);
            f_geo = (nxJ[i,e], nyJ[i,e], sJ[i,e], cf[i,e]);
            FS1[i,e], FS2[i,e], FS3[i,e] = swe_2d_ID_n_surface(ul, ur, du, f_geo, g);
        end
    end
    # # TODO: make faster
    mul!(f1_ID, Pf, FS1); mul!(f2_ID, Pf, FS2); mul!(f3_ID, Pf, FS3);
    # @batch for e = 1:nelem
    #     for i = 1:Nfq # loop over surface nodes
    #         j = Fmask[i] # face node corresponding to volume node
    #         f1_ID[i, e] += wf[j] * FS1[j, e]
    #         f2_ID[i, e] += wf[j] * FS2[j, e]
    #         f3_ID[i, e] += wf[j] * FS3[j, e]
    #     end
    # end

    # @. rhs1_ID = 0;
    fill!(rhs1_ID, zero(eltype(rhs1_ID)));
    #build low order solution first
    @batch for e = 1:nelem
        rxJ_i = rxJ[1,e]; sxJ_i = sxJ[1,e];
        ryJ_i = ryJ[1,e]; syJ_i = syJ[1,e];
        vgeo_e = (rxJ_i, sxJ_i, ryJ_i, syJ_i);
        for i=1:Nq
            UL_E = (h[i,e], hu[i,e], hv[i,e]);
            for j=i:Nq
                UR_E = (h[j,e], hu[j,e], hv[j,e])
                cij = C[i,j,e]
                if i!= j
                    fv1_i_ES, fv1_j_ES = swe_2d_ES_h(UL_E, UR_E, Qr_ES, Qs_ES, vgeo_e, i, j, g);
                    rhs1_ES[i,e] += fv1_i_ES; rhs1_ES[j,e] += fv1_j_ES;
                end
                if cij!=0
                    fv1_i_ID = swe_2d_ID_h(UR_E, Qr_ID, Qs_ID, vgeo_e, i, j, g);
                    fv1_j_ID = swe_2d_ID_h(UL_E, Qr_ID, Qs_ID, vgeo_e, j, i, g);

                    rhs1_ID[i,e] += fv1_i_ID; rhs1_ID[j,e] += fv1_j_ID;
                    lambda_i = abs(u[i,e]*C_x[i,j,e]+v[i,e]*C_y[i,j,e])+sqrt(g*h[i,e])
                    lambda_j = abs(u[j,e]*C_x[j,i,e]+v[j,e]*C_y[j,i,e])+sqrt(g*h[j,e])
                    lambda = max(lambda_i, lambda_j)
                    d1 = co_d * cij * lambda * (h[j,e]  - h[i,e]);
                    rhs1_ID[i,e] -= d1; rhs1_ID[j,e] += d1;
                end
            end
        end
    end
    @. rhs1_ID = rhs1_ID + f1_ID;
    # h_L_next = h - M_inv * rhs1_ID *dt;
    # mul!(Mrhs1_ID, M_inv, rhs1_ID);
    @batch for e = 1:nelem
        J_e = J[1,e];
        for i = 1:Nfq
            Mrhs1_ID[i,e] = M_inv[i,i]*rhs1_ID[i,e]/J_e;
        end
    end
    # axpy!(a, X, Y) Overwrite Y with X*a + Y, where a is a scalar. Return Y.
    lmul!(dt, Mrhs1_ID)
    @. h_L_next = h - Mrhs1_ID;

    @. rhs1_ES = rhs1_ES + f1_ES;
    # mul!(Mrhs1_ES, M_inv, rhs1_ES);
    @batch for e = 1:nelem
        J_e = J[1,e];
        for i = 1:Nfq
            Mrhs1_ES[i,e] = M_inv[i,i]*rhs1_ES[i,e]/J_e;
        end
    end
    lmul!(dt, Mrhs1_ES)
    @. h_H_next = h - Mrhs1_ES;

    mul!(h_L_next_f, E, h_L_next);
    mul!(f1_IDf, E, f1_ID); mul!(f1_ESf, E, f1_ES);
    @. f1_ESf = f1_ESf - f1_IDf
    @batch for e = 1:nelem
        J_e = J[1,e];
        for i = 1:Nfq
            FS1[i,e] = Mf_inv[i,i]*f1_ESf[i,e]/J_e;
        end
    end
    # mul!(FS1, Mf_inv, f1_ESf);
    lmul!(Nq*dt,FS1);
    @. lf = (h_L_next_f -tol)/FS1;
    # mul!(FS1, Mf_inv, f1_ESf);
    # lmul!(Nq*dt,FS1);
    set_val_com!(f1_ESf, lf, tol, 1);
    set_val_com!(h_L_next_f, lf, tol, 0);
    set_val_com!(hf, lf, tol, 0);

    @batch for e = 1:nelem
        for i = 1:Nfq
            lf[i,e] = min(lf[i,e], lf[mapP[i,e]]);
        end
    end
    @. lf = min(1, lf);
    @. lf = max(lf, 0);
    # lf = E'*lf;
    mul!(lff, E', lf);
    @. lff = 0;
    # @. lff = 1;
    @batch for e = 1:nelem
        # if minimum(h_H_next[:,e]) > tol
        if minimum(view(h_H_next, :, e)) > tol
            lff[:,e] .= 1
        end
        for i=1:Nq
            rhs1_CL[i,e] = f1_ID[i,e] + lff[i,e]*(f1_ES[i,e] - f1_ID[i,e]);
            rhs2_CL[i,e] = f2_ID[i,e] + lff[i,e]*(f2_ES[i,e] - f2_ID[i,e]);
            rhs3_CL[i,e] = f3_ID[i,e] + lff[i,e]*(f3_ES[i,e] - f3_ID[i,e]);
        end
    end
    # @show maximum(abs.(f2_ID)), maximum(abs.(f3_ID))
    # @show maximum(abs.(f2_ES)), maximum(abs.(f3_ES))
    # @show maximum(abs.(rhs2_CL)), maximum(abs.(rhs3_CL))
    # @. rhs1_CL = f1_ID + lff*(f1_ES - f1_ID);
    # @. rhs2_CL = f2_ID + lff*(f2_ES - f2_ID);
    # @. rhs3_CL = f3_ID + lff*(f3_ES - f3_ID);


    #volume part
    # loop over all elements
    @batch for e = 1:nelem
        rxJ_i = rxJ[1,e]; sxJ_i = sxJ[1,e];
        ryJ_i = ryJ[1,e]; syJ_i = syJ[1,e];
        vgeo_e = (rxJ_i, sxJ_i, ryJ_i, syJ_i);
        ## fix here, allocation
        # b_e = btm[:,e];
        # h_e_min = minimum(h_H_next[:,e]);
        ###

        h_e_min = minimum(view(h_H_next,:,e));
        J_e = J[1,e];
        for i=1:Nq
            UL_E = (h[i,e], hu[i,e], hv[i,e]);
            for j=i:Nq
                UR_E = (h[j,e], hu[j,e], hv[j,e])
                fv1_i_ES, fv2_i_ES, fv3_i_ES, fv1_j_ES, fv2_j_ES, fv3_j_ES = swe_2d_esdg_vol(UL_E, UR_E, ops, vgeo_e, i, j, view(btm, :,e), g)

                cij = C[i,j]
                fv1_i_ID = 0.0; fv2_i_ID = 0.0; fv3_i_ID = 0.0;
                fv1_j_ID = 0.0; fv2_j_ID = 0.0; fv3_j_ID = 0.0;
                if h_e_min<tol
                # if true
                    if C[i,j]!=0 || i == j
                        fv1_i_ID, fv2_i_ID, fv3_i_ID = swe_2d_ID_vol(UR_E, ops, vgeo_e, i, j, view(btm, :,e), g);
                        fv1_j_ID, fv2_j_ID, fv3_j_ID = swe_2d_ID_vol(UL_E, ops, vgeo_e, j, i, view(btm, :,e), g);

                        # fv1_i_ID, fv2_i_ID, fv3_i_ID = swe_2d_ID_vol(UR_E, Qr_ID, Qs_ID, vgeo_e, i, j);
                        # fv1_j_ID, fv2_j_ID, fv3_j_ID = swe_2d_ID_vol(UL_E, Qr_ID, Qs_ID, vgeo_e, j, i);
                        lambda_i = abs(u[i,e]*C_x[i,j,e]+v[i,e]*C_y[i,j,e])+sqrt(g*h[i,e])
                        lambda_j = abs(u[j,e]*C_x[j,i,e]+v[j,e]*C_y[j,i,e])+sqrt(g*h[j,e])
                        lambda = max(lambda_i, lambda_j)
                        # d1 = 0; d2 = 0; d3 = 0
                        # if h[i,e]>tol &&  h[i,e]>tol
                        # d1 = cij * lambda * (h[j,e] + b_e[j]  - h[i,e] - b_e[i]);
                        d1 = co_d*cij * lambda * (h[j,e] - h[i,e]);
                        # if h[i,e]<=tol ||  h[j,e]<=tol
                        #     d1 = 0;
                        # end
                        # d1 = 0;
                        d2 = co_d*cij * lambda * (hu[j,e] - hu[i,e]);
                        d3 = co_d*cij * lambda * (hv[j,e] - hv[i,e]);
                        # end
                        fv1_i_ID -= d1
                        fv2_i_ID -= d2
                        fv3_i_ID -= d3
                        fv1_j_ID += d1
                        fv2_j_ID += d2
                        fv3_j_ID += d3
                    end

                    l_ij = (h_L_next[i,e] -tol)/(Nq*M_inv[i,i]/J_e*(fv1_i_ES-fv1_i_ID)*dt);
                    # l_ij = (h_L_next[i,e] -tol)/(wq[j]/2*M_inv[i,i]*(fv1_i_ES-fv1_i_ID)*dt);
                    if fv1_i_ES-fv1_i_ID<tol
                        l_ij = 1.0
                    end
                    l_ji = (h_L_next[j,e] -tol)/(Nq*M_inv[j,j]/J_e*(fv1_j_ES-fv1_j_ID)*dt);
                    # l_ji = (h_L_next[j,e] -tol)/(wq[i]/2*M_inv[j,j]*(fv1_j_ES-fv1_j_ID)*dt);
                    if fv1_j_ES-fv1_j_ID<tol
                        l_ji = 1.0
                    end
                    l = min(l_ij, l_ji);
                    l = min(1.0,l);
                    l = max(l,0.0);

                    if (h[i,e] < tol) || (h[j,e] < tol) || (h_L_next[i,e]< tol) || (h_L_next[j,e]< tol) #|| fv1_ES-fv1_i_ID < tol
                        l = 0.0;
                    end

                    # if l <1 && h[i,e]>4
                    #     @show h[i,e], i, e,l
                    #     @show h_L_next[i,e], fv1_i_ES, fv1_i_ID
                    #     @show h_L_next[j,e], fv1_j_ES, fv1_j_ID
                    #     error("limiting in force")
                    # end
                else
                    l = 1;
                end
                # l = 1.0;
                l = 0.0;
                rhs1_CL[i,e] += fv1_i_ID + l * (fv1_i_ES-fv1_i_ID);
                rhs2_CL[i,e] += fv2_i_ID + l * (fv2_i_ES-fv2_i_ID);
                rhs3_CL[i,e] += fv3_i_ID + l * (fv3_i_ES-fv3_i_ID);
                if i!= j
                    rhs1_CL[j,e] += fv1_j_ID + l * (fv1_j_ES-fv1_j_ID);
                    rhs2_CL[j,e] += fv2_j_ID + l * (fv2_j_ES-fv2_j_ID);
                    rhs3_CL[j,e] += fv3_j_ID + l * (fv3_j_ES-fv3_j_ID);
                end
            end
        end
    end
    # @show maximum(abs.(rhs2_CL)), maximum(abs.(rhs3_CL))
    @batch for e = 1:nelem
        J_e = J[1,e];
        for i = 1:Nq
            rhs1_CL[i,e] = -M_inv[i,i]*rhs1_CL[i,e]/J_e;
            rhs2_CL[i,e] = -M_inv[i,i]*rhs2_CL[i,e]/J_e;
            rhs3_CL[i,e] = -M_inv[i,i]*rhs3_CL[i,e]/J_e;
        end
    end
    # @show maximum(abs.(rhs2_CL)), maximum(abs.(rhs3_CL))
    return rhs1_CL, rhs2_CL, rhs3_CL
end