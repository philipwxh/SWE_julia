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

# g = 9.8
# "Approximation parameters"
# N   = 1 # The order of approximation
# CFL = 1/16
# T   = 60 # endtime
# MAXIT = 100#00000
# @show N, MAXIT, T
# ts_ft= 4
# tol = 1e-3
# global n_mr = 0.03
# global gm_bf = 4/3
# # t_plot = LinRange(0, 1,21);
# t_plot = [T]
# save_data = false
# plot_all = false
# plot_last = true
# plot_init = false

# rd = RefElemData(Tri(), SBP{Kubatko{LegendreFaceNodes}}(), N);
# (VXY), EToV = readGmsh2D("malpasset.msh");
# md = MeshData(VXY, EToV, rd);

# # # Construct matrices on reference elements
# @unpack r,s,rf,sf,wf,rq,sq,wq,nrJ,nsJ = rd
# @unpack VDM,V1,Vq,Vf,Dr,Ds,M,Pq,LIFT = rd

# Qr_ES = M*Dr;
# Qs_ES = M*Ds;
# Pq = M\(Vq'*diagm(wq));

# # diff. matrices redefined in terms of quadrature points
# Qr_ES = Pq'*Qr_ES*Pq;
# Qs_ES = Pq'*Qs_ES*Pq;
# E_ES = Vf*Pq;

# # # Need to choose α so that Qr, Qs have zero row sums (and maybe a minimum number of neighbors)
# # # α = 4 # for N=1
# # # α = 2.5 #for N=2
# α = 4 # for N=3
# Qr_ID,Qs_ID,E,Br,Bs,A = build_meshfree_sbp(rq,sq,wq,rf,sf,wf,nrJ,nsJ,α)
# E = floor.(E.+0.1);
# if (norm(sum(Qr_ID,dims=2)) > 1e-10) | (norm(sum(Qs_ID,dims=2)) > 1e-10)
#     error("Qr_ID or Qs_ID doesn't sum to zero for α = $α")
# end
# Qr_ID = Matrix(droptol!(sparse(Qr_ID),1e-15))
# Qs_ID = Matrix(droptol!(sparse(Qs_ID),1e-15))
# Qrskew_ID = .5*(Qr_ID-transpose(Qr_ID))
# Qsskew_ID = .5*(Qs_ID-transpose(Qs_ID))

# @unpack x, y, xf, yf, mapP, mapM, mapB, num_elements = md

# @unpack rxJ, sxJ, ryJ, syJ, J, sJ, nxJ, nyJ = md
# nx = nxJ./sJ; ny = nyJ./sJ;
# QNr = [Qr_ES - .5*E_ES'*Br*E_ES .5*E_ES'*Br;
#     -.5*Br*E_ES .5*Br];
# QNs = [Qs_ES - .5*E_ES'*Bs*E_ES .5*E_ES'*Bs;
#     -.5*Bs*E_ES .5*Bs];

# VN_sbp = [Matrix{Float64}(I, length(wq), length(wq)); E]
# QNr_sbp = VN_sbp'*QNr*VN_sbp;
# # # @show norm(QNr_sbp+QNr_sbp' - E'*diagm(wf.*nrJ)*E)
# QNs_sbp = VN_sbp'*QNs*VN_sbp;
# # # @show norm(QNs_sbp+QNs_sbp' - E'*diagm(wf.*nsJ)*E)

# Qrskew_ES = .5*(QNr_sbp-QNr_sbp');
# Qsskew_ES = .5*(QNs_sbp-QNs_sbp');

# M_inv = diagm(@. 1/(wq))
# M_inv_neg = -M_inv;
# Mf_inv = E*M_inv*E';
# Pf = transpose(E)*diagm(wf)
# Cf = abs.(diag(Qr_ID)[1:length(wf)]*transpose( rxJ[1,:] + ryJ[1,:] )#diag(Qr_ID)[1:length(wf)] 
#         + diag(Qs_ID)[1:length(wf)]*transpose( sxJ[1,:] + syJ[1,:] ))#*diag(Qs_ID)[1:length(wf)])



# cij_x = Array{Float64}(undef,size(Qs_ID, 1),size(Qs_ID, 1),num_elements);
# cij_y = Array{Float64}(undef,size(Qs_ID, 1),size(Qs_ID, 1),num_elements);
# for i = 1:num_elements
#     cij_x[:,:, i] =  rxJ[1, i]*Qr_ID + sxJ[1,i]*Qs_ID
#     cij_y[:,:, i] =  ryJ[1, i]*Qr_ID + syJ[1,i]*Qs_ID
# end
# C = sqrt.(cij_x.*cij_x+cij_y.*cij_y)
# C_x = cij_x./C; C_y = cij_y./C
# replace!(C_x, NaN=>0); replace!(C_y, NaN=>0);
# J_max = maximum(abs.(J))
# J_min = minimum(abs.(J))
# "Time integration"
# rk4a,rk4b,rk4c = ck45()
# CN = (N+1)*(N+2)/2  # estimated trace constant

# # dT = CFL * 2 / (CN*num_elements)
# dT = CFL * minimum(abs.(J[1:size(sJ,1), :]./sJ))/ CN * 2
# Nsteps = convert(Int,ceil(T/dT))
# dt = T/Nsteps

# # "initial conditions"
# xq = Vq*x
# yq = Vq*y
# btm_data = readdlm("malpasset_dam_bathymetry.txt", Float64, skipstart=1);
# itp = interpolate(Multiquadratic(), btm_data[:,1:2]', btm_data[:,3]);
# itpNN = interpolate(NearestNeighbor(), btm_data[:,1:2]', btm_data[:,3]);
# btm = evaluate(itpNN, [md.VX  md.VY]')
# inp_quad = (vandermonde(Tri(), 1, nodes(Tri(), N)...) / vandermonde(Tri(), 1, nodes(Tri(), 1)...));
# btm_q = rd.V1*btm[EToV']

# h = xq*0 .+tol;
# dam_x = [4701.18,4656.5]; 
# dam_y = [4143.41,4392.1];
# k = (dam_y[2]-dam_y[1])/(dam_x[2]-dam_x[1]);
# flag = yq .- dam_y[1] - k*( xq .- dam_x[1] );
# water_idx = findall(x->x<=0, flag);
# h[water_idx] = 100 .- btm_q[water_idx] .+ tol;

# bnd_y_idx = findall(y->abs( y .- 5250 ) < 150, yq);
# bnd_x_idx = findall(x->abs( x .- 4500 ) < 200, xq);
# bnd_idx = bnd_x_idx[findall(x->x in bnd_y_idx, bnd_x_idx)]
# h[bnd_idx] .= tol;

# ocean_idx = findall(x->x<=0, btm_q);

# h[ocean_idx] = tol .- btm_q[ocean_idx];

# # btm_q = btm_q*0;
# #test wb
# # h = 110.0 .-btm_q;
# # @show minimum(h)
# if plot_init
#     gr(aspect_ratio=1,legend=false,
#     markerstrokewidth=0,markersize=2)

#     # "plotting nodes"
#     @unpack rp,sp, Vp = rd
#     vv = Vp*Pq*h

#     fig1 = Makie.Figure();
#     fig2 = Makie.Figure();
#     ax1 = Makie.Axis(fig1[1, 1],
#                     aspect = DataAspect(),
#                     show_axis=false, resolution = (2500,2500));
#     ax2 = Makie.Axis3(fig2[1,1], aspect = (1, 1, 1),
#                     elevation = .25*pi, azimuth = -.25*pi,
#                     show_axis=false, resolution = (2500, 2500));
#     plot_data = Vp*Pq*(h+btm_q)

#     plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y), set_z_coordinate_to_zero=true),
#             color=vec(plot_data),
#             # shading = false,
#             colormap =:blues,
#             # colormap =:viridis,
#             );
#     # plt = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
#     #         color=vec(plot_data), shading = false, colormap = :blues);
#     #
#     # ax = [ax1, ax2]
#     Makie.hidespines!(ax1)
#     Makie.hidedecorations!(ax1)
#     # save("h1.0_mf_90.png", fig, px_per_unit = 2)
#     Makie.tightlimits!(ax1)
#     Makie.Colorbar(fig1[1,2], plt1)
#     # trim!(fig.layout);
#     rowsize!(fig1.layout, 1, ax1.scene.px_area[].widths[2]);
#     # rowsize!(fig.layout[1,3], 1, ax1.scene.px_area[].widths[2]*1.2)
#     # rowsize!(fig.layout, 1, Aspect(3, 1))
#     fig1

#     plt2 = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
#         color=vec(plot_data),
#         # shading = false,
#         colormap = :blues
#         # colormap = :viridis
#         );
#     Makie.hidespines!(ax2)
#     Makie.hidedecorations!(ax2)
#     fig2
# end
# h0 = copy(h)
# hu = h*0.0;
# hv = h*0.0;
# # error("end here")
# "pack arguments into tuples"
# Fmask = [findfirst(@. abs(rf[i] - rq) + abs(sf[i] - sq) < 100*eps()) for i in eachindex(rf)]
# Nq = size(h,1); Nfq = size(E,1);nelem = size(h,2);
# ops = ( Qrskew_ID,Qsskew_ID, Qr_ID, Qs_ID,
#         Qrskew_ES,Qsskew_ES, QNr_sbp, QNs_sbp,
#         E, M, M_inv, Mf_inv, Pf, Fmask);
# dis_cst = (Cf, C, C_x, C_y, Nq, Nfq, nelem)
# vgeo = (rxJ,sxJ,ryJ,syJ,J)
# fgeo = (nxJ,nyJ,sJ, nx, ny)
# nodemaps = (mapP,mapB)
# (gQNxb, gQNyb) =  ESDG_bottom(QNr_sbp, QNs_sbp, btm_q, vgeo, g);
# # U = (h, hu, hv, btm, gQNxb, gQNyb)
# U = (h, hu, hv, btm_q)

# rhs1_ID = zeros(Float64, size(h)); Mrhs1_ID = zeros(Float64, size(h));
# rhs1_ES = zeros(Float64, size(h)); Mrhs1_ES = zeros(Float64, size(h));
# lf = zeros(Float64, size(E*h));lff = zeros(Float64, size(h));
# f1_ES = zeros(Float64, size(h)); f2_ES = zeros(Float64, size(h)); f3_ES= zeros(Float64, size(h));
# f1_ID = zeros(Float64, size(h)); f2_ID = zeros(Float64, size(h)); f3_ID= zeros(Float64, size(h));
# rhs1_CL = zeros(Float64, size(h)); rhs2_CL = zeros(Float64, size(h)); rhs3_CL = zeros(Float64, size(h));
# u = zeros(Float64, size(h)); v = zeros(Float64, size(h));uf = zeros(Float64, size(E*h)); vf = zeros(Float64, size(E*h));
# hf = zeros(Float64, size(E*h)); huf = zeros(Float64, size(E*h));hvf = zeros(Float64, size(E*h));
# hP = zeros(Float64, size(E*h)); huP = zeros(Float64, size(E*h));hvP = zeros(Float64, size(E*h));
# dh = zeros(Float64, size(E*h)); dhu = zeros(Float64, size(E*h));dhv = zeros(Float64, size(E*h));
# lambdaf = zeros(Float64, size(E*h)); lambdaP = zeros(Float64, size(E*h)); cf = zeros(Float64, size(E*h));
# h_L_next = zeros(Float64, size(h)); h_L_next_f = zeros(Float64, size(E*h))
# h_H_next = zeros(Float64, size(h)); h_H_next_f = zeros(Float64, size(E*h))
# f1_IDf = zeros(Float64, size(E*h)); f1_ESf = zeros(Float64, size(E*h));
# rhs_1 = (zeros(Float64, size(h)), zeros(Float64, size(h)), zeros(Float64, size(h)))
# rhs_2 = (zeros(Float64, size(h)), zeros(Float64, size(h)), zeros(Float64, size(h)))
# htmp = zeros(Float64, size(h)); hutmp = zeros(Float64, size(h)); hvtmp = zeros(Float64, size(h));
# FS1 = zeros(Float64, size(E*h)); FS2 = zeros(Float64, size(E*h)); FS3 = zeros(Float64, size(E*h));
# UL_E = (Float64, Float64, Float64); UR_E = (Float64, Float64, Float64); vgeo_e =(Float64, Float64, Float64, Float64);
# pre_allo = (u, v, uf, vf, hf, huf, hvf, hP, huP, hvP, dh, dhu, dhv,
#             rhs1_ID, Mrhs1_ID, rhs1_ES, Mrhs1_ES, lf, lff, f1_ES, f2_ES, f3_ES, f1_ID, f2_ID, f3_ID,
#             rhs1_CL, rhs2_CL, rhs3_CL, lambdaf, lambdaP, cf, UL_E, UR_E, vgeo_e,
#             h_L_next, h_L_next_f, f1_IDf, f1_ESf, FS1, FS2, FS3)


function swe_2d_rhs(U,ops,dis_cst,vgeo,fgeo,nodemaps, dt, tol, g, pre_allo)
    # unpack args
    h, hu, hv, btm = U
    (Qr_ID,Qs_ID, Qrb_ID, Qsb_ID,Qr_ES,Qs_ES,QNr_sbp, QNs_sbp,
     E, M, M_inv, Mf_inv, Pf,Fmask ) = ops
    Cf, C, C_x, C_y, Nq, Nfq, nelem= dis_cst
    rxJ,sxJ,ryJ,syJ,J = vgeo
    nxJ,nyJ,sJ,nx,ny = fgeo
    mapP,mapB = nodemaps
    (u, v, uf, vf, hf, huf, hvf, hP, huP, hvP, dh, dhu, dhv,
    rhs1_ID, Mrhs1_ID, rhs1_ES, Mrhs1_ES, lf, lff,
    f1_ES, f2_ES, f3_ES, f1_ID, f2_ID, f3_ID, rhs1_CL, rhs2_CL, rhs3_CL,
    lambdaf, lambdaP, cf, UL_E, UR_E, vgeo_e,
    h_L_next, h_L_next_f, f1_IDf, f1_ESf, FS1, FS2, FS3)= pre_allo
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
        vu_dot_n  = huP[i]*nx[i]+ hvP[i]*ny[i];
        # vu_dot_n = vu_dot_n;
        huP[i] =  huP[i] - 2*vu_dot_n*nx[i];
        hvP[i] =  hvP[i] - 2*vu_dot_n*ny[i];
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
                    d1 = cij * lambda * (h[j,e]  - h[i,e]);
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
    # @. lff = 0;
    @. lff = 1;
    @batch for e = 1:nelem
        if minimum(h_H_next[:,e]) > tol
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
        b_e = btm[:,e];
        h_e_min = minimum(h_H_next[:,e]);
        J_e = J[1,e];
        for i=1:Nq
            UL_E = (h[i,e], hu[i,e], hv[i,e]);
            for j=i:Nq
                UR_E = (h[j,e], hu[j,e], hv[j,e])
                fv1_i_ES, fv2_i_ES, fv3_i_ES, fv1_j_ES, fv2_j_ES, fv3_j_ES = swe_2d_esdg_vol(UL_E, UR_E, ops, vgeo_e, i, j, b_e, g)

                cij = C[i,j,e]
                fv1_i_ID = 0.0; fv2_i_ID = 0.0; fv3_i_ID = 0.0;
                fv1_j_ID = 0.0; fv2_j_ID = 0.0; fv3_j_ID = 0.0;
                # if h_e_min<tol
                if true
                    if C[i,j]!=0 || i == j
                        fv1_i_ID, fv2_i_ID, fv3_i_ID = swe_2d_ID_vol(UR_E, ops, vgeo_e, i, j, b_e, g);
                        fv1_j_ID, fv2_j_ID, fv3_j_ID = swe_2d_ID_vol(UL_E, ops, vgeo_e, j, i, b_e, g);

                        # fv1_i_ID, fv2_i_ID, fv3_i_ID = swe_2d_ID_vol(UR_E, Qr_ID, Qs_ID, vgeo_e, i, j);
                        # fv1_j_ID, fv2_j_ID, fv3_j_ID = swe_2d_ID_vol(UL_E, Qr_ID, Qs_ID, vgeo_e, j, i);
                        lambda_i = abs(u[i,e]*C_x[i,j,e]+v[i,e]*C_y[i,j,e])+sqrt(g*h[i,e])
                        lambda_j = abs(u[j,e]*C_x[j,i,e]+v[j,e]*C_y[j,i,e])+sqrt(g*h[j,e])
                        lambda = max(lambda_i, lambda_j)
                        # d1 = 0; d2 = 0; d3 = 0
                        # if h[i,e]>tol &&  h[i,e]>tol
                        # d1 = cij * lambda * (h[j,e] + b_e[j]  - h[i,e] - b_e[i]);
                        d1 = cij * lambda * (h[j,e] - h[i,e]);
                        # if h[i,e]<=tol ||  h[j,e]<=tol
                        #     d1 = 0;
                        # end
                        # d1 = 0;
                        d2 = cij * lambda * (hu[j,e] - hu[i,e]);
                        d3 = cij * lambda * (hv[j,e] - hv[i,e]);
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
            v_l2 = sqrt( (hu[i,e]/h[i,e])^2 + (hv[i,e]/h[i,e])^2 );
            h_modified = max(h[i,e]^gm_bf, 2*g*n_mr^2*dt*v_l2);
            rhs2_CL[i,e] += (2*g*n_mr^2*hu[i,e]*v_l2*M[i,i])/(h[i,e]^gm_bf+h_modified)
            rhs3_CL[i,e] += (2*g*n_mr^2*hv[i,e]*v_l2*M[i,i])/(h[i,e]^gm_bf+h_modified)
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

DT = zeros(MAXIT)
t = 0.0
pl_idx = 1;
global i;
# filename = string("h",string(N),"_",string(K1D),"_","0_dam_break.csv");
# writedlm( filename,  h, ',');
# allocate = @allocated begin

@time begin
for i = 1:MAXIT
    # if i%1000 == 0
    if i%10 == 0
        @show maximum(h), maximum(abs.(u)), maximum(abs.(v))
        @show i, t
    end
    # @show i, t
    global h, hu, hv, U, t, t_plot, pl_idx, dt, dt1, dt2, lambda, filename, J_max, J_min
    global U,t, t_plot, pl_idx, dt, filename, htmp, hutmp, hvtmp, utmp, rhs_1, rhs_2, h_max
    # Heun's method - this is an example of a 2nd order SSP RK method
    # local rhs_ES1, rhs_ID1 = swe_2d_rhs(u,ops,dis_cst,vgeo,fgeo,nodemaps, dt1)
    # lambda = maximum(sqrt.((hu./h).^2+(hv./h).^2)+sqrt.(g.*h))
    @. u = hu/h; @. v = hv/h;
    @. u = sqrt((u)^2+(v)^2)+sqrt(g*h)
    lambda = maximum(u)
    dt1 = min(min(T,t_plot[pl_idx])-t, minimum(wq)*J_min/(ts_ft*lambda), dT);
    # dt1 = dT
    rhs_1 = swe_2d_rhs(U,ops,dis_cst,vgeo,fgeo,nodemaps, dt1, tol, g, pre_allo)
    # @show lambda, dt1
    
    @. htmp  = h  + dt1*rhs_1[1]
    @. hutmp = hu + dt1*rhs_1[2]
    @. hvtmp = hv + dt1*rhs_1[3]
    h_max = maximum(htmp)
    # utmp = (htmp, hutmp)
    zero_vel!(htmp, hutmp, hvtmp, tol);
    # vel_reg!(htmp, hutmp, hvtmp, 55, tol);
    # view(hutmp, findall(x->x<2*tol, htmp)) .= 0

    h_min, pos = findmin(htmp)
    if h_min <= 0
       # @show L1
       error("htmp_min<0 ", h_min, pos, "iteration ", i )
    end
    @. u = hutmp/htmp; @. v = hvtmp/htmp;
    # @show maximum(abs.(u)), maximum(abs.(v))
    @. u = sqrt((u)^2+(v)^2)+sqrt(g*htmp)
    lambda = maximum(u)
    # @show lambda, dt1
    dt2 = min(min(T,t_plot[pl_idx])-t, minimum(wq)*J_min/(ts_ft*lambda), dT);
    # dt2 = dT
    # @show dt1, dt2
    while dt2<dt1
        dt1 = dt1/2
        @. htmp  = h  + dt1*rhs_1[1]
        @. hutmp = hu + dt1*rhs_1[2]
        @. hvtmp = hv + dt1*rhs_1[3]
        # utmp = (htmp, hutmp)
        zero_vel!(htmp, hutmp, hvtmp, tol);
        # vel_reg!(htmp, hutmp, hvtmp, 55, tol);
        @. u = hutmp/htmp; @. v = hv/htmp;
        @. u = sqrt((u)^2+(v)^2)+sqrt(g*h)
        lambda = maximum(u)
        dt2 = min(min(T,t_plot[pl_idx])-t, minimum(wq)*J_min/(ts_ft*lambda), dT);
        # @show dt2
    end
    # @show dt1, dt2
    # zero_vel!(htmp, hutmp, hvtmp, tol);
    utmp = (htmp, hutmp, hvtmp, btm_q)
    dt = min(dt1, dt2)
    rhs_2 = swe_2d_rhs(utmp,ops,dis_cst,vgeo,fgeo,nodemaps, dt, tol, g, pre_allo)
    # dt = min(dt1, dt2)
    # @show i, dt, t
    # s1 = sum(h)+sum(hu)+sum(hv)
    @. h  += .5*dt*(rhs_1[1] + rhs_2[1])
    @. hu += .5*dt*(rhs_1[2] + rhs_2[2])
    @. hv += .5*dt*(rhs_1[3] + rhs_2[3])
    zero_vel!(h, hu, hv, tol);
    # vel_reg!(h, hu, hv, 55, tol);
    
    # @show L2, h, dt
    h_min, pos = findmin(h)

    if h_min <=0
        # @show L2
        @show maximum(hu./h)
        error("h_min<0 ", h_min, pos, "iteration ", i )
    end
    U = (h,hu,hv, btm_q)
    t += dt
    # @show L1, L2
    DT[i] = dt
    i +=1
    if t>= t_plot[pl_idx] #|| i==Nsteps
        # @show t
        if save_data
            filename = string("malpasset_h",string(N),"_",pl_idx,".csv");
            writedlm( filename,  h, ',');
            filename = string("malpasset_hu",string(N),"_",string(K1D),"_",pl_idx,".csv");
            writedlm( filename,  hu, ',');
            filename = string("malpasset_hv",string(N),"_",string(K1D),"_",pl_idx,".csv");
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
    vv = Vp*Pq*h0

    fig1 = Makie.Figure();
    fig2 = Makie.Figure();
    ax1 = Makie.Axis(fig1[1, 1],
                    aspect = DataAspect(),
                    show_axis=false, resolution = (2500,2500));
    ax2 = Makie.Axis3(fig2[1,1], aspect = (1, 1, 1),
                    elevation = .25*pi, azimuth = -.25*pi,
                    show_axis=false, resolution = (2500, 2500));
    plot_data = Vp*Pq*(h0)

    plt1 = Makie.mesh!(ax1, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y), set_z_coordinate_to_zero=true),
            color=vec(plot_data),
            # shading = false,
            colormap =:blues,
            # colormap =:viridis,
            );
    xb = xf[mapB]; yb = yf[mapB];
    xh = xq[findall(x->x>80, abs.(h))]; yh = yq[findall(x->x>80, abs.(h))];
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
    for idx in Integer.(LinRange(5, 17,4))
    # idx = 2
    tag = "l_n_H";
    filename = string("db_", tag, "_h3_", string(K1D), "_", string(idx), ".csv");
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
    filename = string("db_", tag, "_h_", string(K1D),"_", string(idx-1), "_90.png");
    save(filename, fig1, px_per_unit = 2)
    # # filename = string("dam_break_low_alpha_h", string(idx*4+1),"_32_45.png");
    filename = string("db_", tag, "_h_", string(K1D),"_", string(idx-1), "_45.png");
    save(filename, fig2, px_per_unit = 2)
    end
end

