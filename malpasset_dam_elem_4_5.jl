using Revise
using Plots
using LinearAlgebra
using SparseArrays
using StaticArrays
using UnPack
using StartUpDG
using DelimitedFiles
using BenchmarkTools
using Polyester
using ScatteredInterpolation
# using Interpolations
include("dg2d_swe_flux.jl")
include("dg2d_swe_mesh_opt.jl")
include("plot_Makie.jl")
include("malpasset_setup.jl")

function swe_2d_rhs(U,ops,dis_cst,vgeo,fgeo,nodemaps, dt, tol, g, pre_allo)
    # unpack args
    h, hu, hv, btm = U
    (Qr_ID,Qs_ID,Qrb_ID,Qsb_ID,Qr_ES,Qs_ES,QNr_sbp, QNs_sbp,
     E, M, M_inv, Mf_inv, Pf,Fmask) = ops
    Cf, C, C_x, C_y, Nq, Nfq, nelem= dis_cst
    rxJ,sxJ,ryJ,syJ,J = vgeo
    nxJ,nyJ,sJ,nx,ny = fgeo
    mapP,mapB = nodemaps
    (u, v, uf, vf, hf, huf, hvf, hP, huP, hvP, dh, dhu, dhv,
    rhs1_ID, Mrhs1_ID, rhs1_ES, Mrhs1_ES, lf, lff,
    f1_ES, f2_ES, f3_ES, f1_ID, f2_ID, f3_ID, rhs1_CL, rhs2_CL, rhs3_CL,
    lambdaf, lambdaP, cf, UL_E, UR_E, vgeo_e,
    h_L_next, h_L_next_f, f1_IDf, f1_ESf, FS1, FS2, FS3)= pre_allo
    co_d = 1
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
            hP[i,e]  = hf[mapP[i,e]];
            huP[i,e] = huf[mapP[i,e]];
            hvP[i,e] = hvf[mapP[i,e]];
        end
    end
    @batch for i in mapB
        vu_dot_n = huf[i]*nx[i] + hvf[i]*ny[i];
        vu_dot_n = vu_dot_n;
        # huP[i] = huf[i] - 2*vu_dot_n*nx[i];
        # hvP[i] = hvf[i] - 2*vu_dot_n*ny[i];
        huP[i] =  -huf[i];
        hvP[i] =  -hvf[i];
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

    mul!(f1_ES, Pf, FS1);mul!(f2_ES, Pf, FS2); mul!(f3_ES, Pf, FS3);

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
                    d1 = co_d*cij * lambda * (h[j,e] - h[i,e]);
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
        for i = 1:Nq
            Mrhs1_ID[i,e] = M_inv[i,i]*rhs1_ID[i,e]/J_e;
        end
    end
    # axpy!(a, X, Y) Overwrite Y with X*a + Y, where a is a scalar. Return Y.
    lmul!(dt, Mrhs1_ID)
    @. h_L_next = h - Mrhs1_ID;
    while  minimum(h_L_next) < tol
        dt = dt/2;
        lmul!(dt, Mrhs1_ID)
        @. h_L_next = h - Mrhs1_ID;
        # @show minimum(h_L_next), dt
        if dt < tol
            @show minimum(h_L_next), dt
            error("dt is too small")
        end
    end

    @. rhs1_ES = rhs1_ES + f1_ES;
    # mul!(Mrhs1_ES, M_inv, rhs1_ES);
    @batch for e = 1:nelem
        J_e = J[1,e];
        for i = 1:Nq
            Mrhs1_ES[i,e] = M_inv[i,i]*rhs1_ES[i,e]/J_e;
        end
    end
    lmul!(dt, Mrhs1_ES)
    @. h_H_next = h - Mrhs1_ES;
   
    @batch for e = 1:nelem
        for i=1:Nq
            rhs1_CL[i,e] = f1_ID[i,e] + f1_ES[i,e];
            rhs2_CL[i,e] = f2_ID[i,e] + f2_ES[i,e];
            rhs3_CL[i,e] = f3_ID[i,e] + f3_ES[i,e];
        end
    end
    
    #volume part
    # loop over all elements
    @batch for e = 1:nelem
        rxJ_i = rxJ[1,e]; sxJ_i = sxJ[1,e];
        ryJ_i = ryJ[1,e]; syJ_i = syJ[1,e];
        vgeo_e = (rxJ_i, sxJ_i, ryJ_i, syJ_i);
        l_e = 1.0;
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

                cij = C[i,j,e]
                fv1_i_ID = 0.0; fv2_i_ID = 0.0; fv3_i_ID = 0.0;
                fv1_j_ID = 0.0; fv2_j_ID = 0.0; fv3_j_ID = 0.0;
                if h_e_min<tol
                # if true
                # if false
                    if cij!=0 || i == j
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
                        d1 = co_d * cij * lambda * (h[j,e] - h[i,e]);
                        # if h[i,e]<=tol ||  h[j,e]<=tol
                        #     d1 = 0;
                        # end
                        # d1 = 0;
                        d2 = co_d * cij * lambda * (hu[j,e] - hu[i,e]);
                        d3 = co_d * cij * lambda * (hv[j,e] - hv[i,e]);
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
                    l = 1.0;
                end

                l_e = min(l,l_e)
                # l_e = 0
                if l_e < tol
                    l_e = 0.0;
                    break
                end
            end
            if l_e < tol
                l_e = 0.0;
                break
            end
        end
        for i=1:Nq
            UL_E = (h[i,e], hu[i,e], hv[i,e]);
            for j=i:Nq
                UR_E = (h[j,e], hu[j,e], hv[j,e])
                fv1_i_ES, fv2_i_ES, fv3_i_ES, fv1_j_ES, fv2_j_ES, fv3_j_ES = swe_2d_esdg_vol(UL_E, UR_E, ops, vgeo_e, i, j, view(btm, :,e), g)

                cij = C[i,j,e]
                fv1_i_ID = 0.0; fv2_i_ID = 0.0; fv3_i_ID = 0.0;
                fv1_j_ID = 0.0; fv2_j_ID = 0.0; fv3_j_ID = 0.0;
                if h_e_min<tol
                # if true
                    if cij!=0 || i == j
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
                        d1 = co_d * cij * lambda * (h[j,e] - h[i,e]);
                        # if h[i,e]<=tol ||  h[j,e]<=tol
                        #     d1 = 0;
                        # end
                        # d1 = 0;
                        d2 = co_d * cij * lambda * (hu[j,e] - hu[i,e]);
                        d3 = co_d * cij * lambda * (hv[j,e] - hv[i,e]);
                        # end
                        fv1_i_ID -= d1
                        fv2_i_ID -= d2
                        fv3_i_ID -= d3
                        fv1_j_ID += d1
                        fv2_j_ID += d2
                        fv3_j_ID += d3
                    end
                end
                # l_e = 1.0;
                # l_e = 0.0;
                rhs1_CL[i,e] += fv1_i_ID + l_e * (fv1_i_ES-fv1_i_ID);
                rhs2_CL[i,e] += fv2_i_ID + l_e * (fv2_i_ES-fv2_i_ID);
                rhs3_CL[i,e] += fv3_i_ID + l_e * (fv3_i_ES-fv3_i_ID);
                if i!= j
                    rhs1_CL[j,e] += fv1_j_ID + l_e * (fv1_j_ES-fv1_j_ID);
                    rhs2_CL[j,e] += fv2_j_ID + l_e * (fv2_j_ES-fv2_j_ID);
                    rhs3_CL[j,e] += fv3_j_ID + l_e * (fv3_j_ES-fv3_j_ID);
                end
            end
            v_l2 = sqrt( (hu[i,e]/h[i,e])^2 + (hv[i,e]/h[i,e])^2 );
            h_modified = max(h[i,e]^gm_bf, 2*g*n_mr^2*dt*v_l2);
            rhs2_CL[i,e] += (2*g*n_mr^2*hu[i,e]*v_l2*M[i,i])*J_e/(h[i,e]^gm_bf+h_modified)
            rhs3_CL[i,e] += (2*g*n_mr^2*hv[i,e]*v_l2*M[i,i])*J_e/(h[i,e]^gm_bf+h_modified)
        end
    end
    # @show maximum(abs.(rhs2_CL)), maximum(abs.(rhs3_CL))
    @batch for e = 1:nelem
        J_e = J[1,e];
        for i = 1:Nq
            rhs1_CL[i,e] = -M_inv[i,i]*rhs1_CL[i,e]/J_e;
            rhs2_CL[i,e] = -M_inv[i,i]*rhs2_CL[i,e]/J_e;
            rhs3_CL[i,e] = -M_inv[i,i]*rhs3_CL[i,e]/J_e;
            # if abs(rhs3_CL[i,e]) > 10000
            #     @show  i, e, h[i,e], hu[i,e], hv[i,e]
            #     error("blow up")
            # end
        end
    end
    # @show maximum(abs.(rhs2_CL)), maximum(abs.(rhs3_CL))
    return rhs1_CL, rhs2_CL, rhs3_CL, dt
end

DT = zeros(MAXIT)
t = 0.0
pl_idx = 1;
global i;

@time begin
for i = 1:MAXIT
    # if i%1000 == 0
    if i%1 == 0
        @show maximum(h), maximum(abs.(u)), maximum(abs.(v))
        @show i, t
    end
    # @show i, t
    global h, hu, hv, U, t, t_plot, pl_idx, dt, dt1, dt2, lambda, filename, J_max, J_min, hh_max, h_min
    global U,t, t_plot, pl_idx, dt, filename, htmp, hutmp, hvtmp, utmp, rhs_1, rhs_2, h_max
    # Heun's method - this is an example of a 2nd order SSP RK method
    # local rhs_ES1, rhs_ID1 = swe_2d_rhs(u,ops,dis_cst,vgeo,fgeo,nodemaps, dt1)
    # lambda = maximum(sqrt.((hu./h).^2+(hv./h).^2)+sqrt.(g.*h))
    hh_max = maximum(h);
    @. u = hu/h; @. v = hv/h;
    @. u = sqrt((u)^2+(v)^2)+sqrt(g*h)
    lambda = maximum(u)

    dt = min(min(T,t_plot[pl_idx])-t, minimum(wq) * J_min / (ts_ft * lambda), dT);
    rhs_1 = swe_2d_rhs(U,ops,dis_cst,vgeo,fgeo,nodemaps, dt, tol, g, pre_allo)
    dt = rhs_1[4]
    @. htmp = h + dt*(rhs_1[1])
    h_min, pos = findmin(htmp)

    while h_min < tol/2
        dt = dt/2
        # @show dt, i
        @. htmp = h + dt*(rhs_1[1])
        h_min, pos = findmin(htmp)
        if dt < tol
            @show dt, i
            error("dt too small")
        end
    end
    @. h  += dt*(rhs_1[1])
    @. hu += dt*(rhs_1[2])
    @. hv += dt*(rhs_1[3])
    zero_vel!(h, hu, hv, tol, hh_max);

    # vel_reg!(htmp, hutmp, hvtmp, 55, tol);
    # @show L2, h, dt
    h_min, pos = findmin(h)

    if h_min <=0
        # @show L2
        @show maximum(hu./h)
        error("h_min<0 ", h_min, pos, "iteration ", i )
    end
    U = (h,hu,hv, btm_q)
    t += dt
    DT[i] = dt
    i +=1
    if t>= t_plot[pl_idx] #|| i==Nsteps
        @show t
        if save_data
            filename = string("malpasset_h",string(N),"_",string(Integer(round(t))),".csv");
            writedlm( filename,  h, ',');
            filename = string("malpasset_hu",string(N),"_",string(Integer(round(t))),".csv");
            writedlm( filename,  hu, ',');
            filename = string("malpasset_hv",string(N),"_",string(Integer(round(t))),".csv");
            writedlm( filename,  hv, ',');
            println("save at iteration",i);
        end
        pl_idx+=1;
        # t_plot += dT*10;
        # println("Number of time steps $i out of $Nsteps")
    end
    if t>=T
        break
    end
end
end
# end; if allocate > 0 println(allocate) end

# DT = DT[1:findmin(DT)[2]-1];
if plot_last
    gr(aspect_ratio=1,legend=false,
    markerstrokewidth=0,markersize=2)

    # "plotting nodes"
    @unpack rp,sp, Vp = rd
    vv = Vp*Pq*h

    fig1 = Makie.Figure();
    fig2 = Makie.Figure();
    ax1 = Makie.Axis(fig1[1, 1],
                    aspect = DataAspect(),
                    show_axis=false, resolution = (2500,2500));
    ax2 = Makie.Axis3(fig2[1,1], aspect = (1, 1, 1),
                    elevation = .25*pi, azimuth = -.25*pi,
                    show_axis=false, resolution = (2500, 2500));
    plot_data = Vp*Pq*(sqrt.((hu./h).^2 + (hv./h).^2))

    plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y), set_z_coordinate_to_zero=true),
            color=vec(plot_data),
            # shading = false,
            colormap =:blues,
            # colormap =:viridis,
            );
    # xb = xf[mapB]; yb = yf[mapB];
    # B_idx = zeros(size(xf))
    # B_idx[mapB] .= 1
    # B_idx = rd.Vf'*B_idx
    # B_idx = findall(x->x>0, idxx);
    xb = xq[B_idx]; yb = yq[B_idx]
    # xb = xq[findall(x->x>90, btm_q)]; yb = yq[findall(x->x>90, btm_q)]
    xh = xq[findall(x->x>100, abs.(h))]; yh = yq[findall(x->x>100, abs.(h))];
    xu = xq[findall(x->x>50, abs.(hu./h))]; yu = yq[findall(x->x>50, abs.(hu./h))];
    xv = xq[findall(x->x>50, abs.(hv./h))]; yv = yq[findall(x->x>50, abs.(hv./h))];
    plt1 = Makie.scatter!(ax1, xb, yb, color = :black, markersize = 2)
    plt1 = Makie.scatter!(ax1, xh, yh, color = :yellow, markersize = 5)
    plt1 = Makie.scatter!(ax1, xu, yu, color = :red, markersize = 5)
    plt1 = Makie.scatter!(ax1, xv, yv, color = :green, markersize = 5)
    # plt = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
    #         color=vec(plot_data), shading = false, colormap = :blues);
    #
    # ax = [ax1, ax2]
    # Makie.hidespines!(ax1)
    # Makie.hidedecorations!(ax1)
    # save("h1.0_mf_90.png", fig, px_per_unit = 2)
    Makie.tightlimits!(ax1)
    Makie.Colorbar(fig1[1,2], plt1)
    # trim!(fig.layout);
    rowsize!(fig1.layout, 1, ax1.scene.px_area[].widths[2]);
    # rowsize!(fig.layout[1,3], 1, ax1.scene.px_area[].widths[2]*1.2)
    # rowsize!(fig.layout, 1, Aspect(3, 1))
    fig1


    plt2 = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
        color=vec(plot_data),
        # shading = false,
        colormap = :blues
        # colormap = :viridis
        );
    # Makie.hidespines!(ax2)
    # Makie.hidedecorations!(ax2)
    fig2
end

if plot_all
    for idx =1:40
        tag = "malpasset_h";
        filename = string(tag, string(N),"_", string(idx), ".csv");
        h  = readdlm(filename, ',', Float64);
        fig1 = Makie.Figure();
        fig2 = Makie.Figure();
        ax1 = Makie.Axis(fig1[1, 1],
                        aspect = DataAspect(),
                        show_axis=false, resolution = (2500,2500));
        ax2 = Makie.Axis3(fig2[1,1], aspect = (1, 1, 1),
                        elevation = .25*pi, azimuth = -.25*pi,
                        show_axis=false, resolution = (2500, 2500));
        plot_data = Vp*Pq*(h)

        plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y), set_z_coordinate_to_zero=true),
                color=vec(plot_data),
                # shading = false,
                colormap =:blues,
                # colormap =:viridis,
                );
        # plt = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
        #         color=vec(plot_data), shading = false, colormap = :blues);
        #
        # ax = [ax1, ax2]
        Makie.hidespines!(ax1)
        Makie.hidedecorations!(ax1)
        # save("h1.0_mf_90.png", fig, px_per_unit = 2)
        Makie.tightlimits!(ax1)
        Makie.Colorbar(fig1[1,2], plt1)
        # trim!(fig.layout);
        rowsize!(fig1.layout, 1, ax1.scene.px_area[].widths[2]);
        # rowsize!(fig.layout[1,3], 1, ax1.scene.px_area[].widths[2]*1.2)
        # rowsize!(fig.layout, 1, Aspect(3, 1))
        # fig1

        plt2 = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
            color=vec(plot_data),
            # shading = false,
            colormap = :blues
            # colormap = :viridis
            );
        Makie.hidespines!(ax2)
        Makie.hidedecorations!(ax2)
        # fig2

        # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_90.png");
        filename = string( tag,"_", string(idx-1), "_90.png");
        save(filename, fig1, px_per_unit = 2)
        # # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_45.png");
        filename = string(tag,"_", string(idx-1),  "_45.png");
        save(filename, fig2, px_per_unit = 2)
    end
end