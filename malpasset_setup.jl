g = 9.8
"Approximation parameters"
N   = 2 # The order of approximation
CFL = 1/16
T   = 100#0 # endtime
MAXIT = 1000000
ts_ft = 4
tol = 1e-4
global n_mr = 0.04
global gm_bf = 4/3
t_plot = LinRange(100, T, Integer(T/100))
# t_plot = [T]
save_data = true
plot_all = false
plot_last = true
plot_init = false
@show N, MAXIT, T, 1/CFL, ts_ft, n_mr

rd = RefElemData(Tri(), SBP{Kubatko{LegendreFaceNodes}}(), N);
(VXY), EToV = readGmsh2D("malpasset.msh");
md = MeshData(VXY, EToV, rd);

# # Construct matrices on reference elements
@unpack r,s,rf,sf,wf,rq,sq,wq,nrJ,nsJ,Nq,Nfq = rd
@unpack VDM,V1,Vq,Vf,Dr,Ds,M,Pq,LIFT = rd

Qr_ES = M*Dr;
Qs_ES = M*Ds;
Pq = M\(Vq'*diagm(wq));

# diff. matrices redefined in terms of quadrature points
Qr_ES = Pq'*Qr_ES*Pq;
Qs_ES = Pq'*Qs_ES*Pq;
E_ES = Vf*Pq;

# # Need to choose α so that Qr, Qs have zero row sums (and maybe a minimum number of neighbors)
# # α = 4 # for N=1
# # α = 2.5 #for N=2
α = 4 # for N=3
Qr_ID,Qs_ID,E,Br,Bs,A = build_meshfree_sbp(rq,sq,wq,rf,sf,wf,nrJ,nsJ,α)
E = floor.(E.+0.1);
if (norm(sum(Qr_ID,dims=2)) > 1e-10) | (norm(sum(Qs_ID,dims=2)) > 1e-10)
    error("Qr_ID or Qs_ID doesn't sum to zero for α = $α")
end
Qr_ID = Matrix(droptol!(sparse(Qr_ID),1e-15))
Qs_ID = Matrix(droptol!(sparse(Qs_ID),1e-15))
Qrskew_ID = .5*(Qr_ID-transpose(Qr_ID))
Qsskew_ID = .5*(Qs_ID-transpose(Qs_ID))

@unpack x, y, xf, yf, mapP, mapM, mapB, num_elements = md

@unpack rxJ, sxJ, ryJ, syJ, J, sJ, nxJ, nyJ = md
nx = nxJ./sJ; ny = nyJ./sJ;
QNr = [Qr_ES - .5*E_ES'*Br*E_ES .5*E_ES'*Br;
    -.5*Br*E_ES .5*Br];
QNs = [Qs_ES - .5*E_ES'*Bs*E_ES .5*E_ES'*Bs;
    -.5*Bs*E_ES .5*Bs];

VN_sbp = [Matrix{Float64}(I, length(wq), length(wq)); E]
QNr_sbp = VN_sbp'*QNr*VN_sbp;
# # @show norm(QNr_sbp+QNr_sbp' - E'*diagm(wf.*nrJ)*E)
QNs_sbp = VN_sbp'*QNs*VN_sbp;
# # @show norm(QNs_sbp+QNs_sbp' - E'*diagm(wf.*nsJ)*E)

Qrskew_ES = .5*(QNr_sbp-QNr_sbp');
Qsskew_ES = .5*(QNs_sbp-QNs_sbp');

M_inv = diagm(@. 1/(wq))
M_inv_neg = -M_inv;
Mf_inv = E*M_inv*E';
Pf = transpose(E)*diagm(wf)
Cf = abs.(diag(Qr_ID)[1:length(wf)]*transpose( rxJ[1,:] + ryJ[1,:] )#diag(Qr_ID)[1:length(wf)] 
        + diag(Qs_ID)[1:length(wf)]*transpose( sxJ[1,:] + syJ[1,:] ))#*diag(Qs_ID)[1:length(wf)])



cij_x = Array{Float64}(undef,size(Qs_ID, 1),size(Qs_ID, 1),num_elements);
cij_y = Array{Float64}(undef,size(Qs_ID, 1),size(Qs_ID, 1),num_elements);
for i = 1:num_elements
    cij_x[:,:, i] =  rxJ[1, i]*Qr_ID + sxJ[1,i]*Qs_ID
    cij_y[:,:, i] =  ryJ[1, i]*Qr_ID + syJ[1,i]*Qs_ID
end
C = sqrt.(cij_x.*cij_x+cij_y.*cij_y)
C_x = cij_x./C; C_y = cij_y./C
replace!(C_x, NaN=>0); replace!(C_y, NaN=>0);
J_max = maximum(abs.(J))
J_min = minimum(abs.(J))
"Time integration"
rk4a,rk4b,rk4c = ck45()
CN = (N+1)*(N+2)/2  # estimated trace constant

# dT = CFL * 2 / (CN*num_elements)
dT = CFL * minimum(abs.(J[1:size(sJ,1), :]./sJ))/ CN * 2
Nsteps = convert(Int,ceil(T/dT))
dt = T/Nsteps

# "initial conditions"
xq = Vq*x
yq = Vq*y
btm_data = readdlm("malpasset_dam_bathymetry.txt", Float64, skipstart=1);
ocean_idx = readdlm("malpasset_dam_ocean_idx.txt", Int64)[:,1];
itp = interpolate(Multiquadratic(), btm_data[:,1:2]', btm_data[:,3]);
itpNN = interpolate(NearestNeighbor(), btm_data[:,1:2]', btm_data[:,3]);
B_idx = zeros(size(xf))
B_idx[mapB] .= 1
B_idx = rd.Vf'*B_idx
B_idx = findall(x->x>0.1, B_idx);
btm = evaluate(itpNN, [md.VX  md.VY]')
inp_quad = (vandermonde(Tri(), 1, nodes(Tri(), N)...) / vandermonde(Tri(), 1, nodes(Tri(), 1)...));
btm_q = rd.V1*btm[EToV']
v_idx = findall(x->x<10000, xq);
sb_idx = []
for e in ocean_idx
    if maximum(btm_q[:,e]) > 0 && minimum(btm_q[:,e]) < 0 
        append!(sb_idx, e)
        btm_q[:,e] .= min.(btm_q[:,e], -tol);
    end
    if maximum(btm_q[:,e]) < 0
        # append!(sb_b_idx, e)
        # # append!(mapB, (e-1)*Nfq+1:e*Nfq)
        mapP[mapP[(e-1)*Nfq+1:e*Nfq]] .= mapM[mapP[(e-1)*Nfq+1:e*Nfq]]
        mapP[(e-1)*Nfq+1:e*Nfq].= (e-1)*Nfq+1:e*Nfq
    end
end

h = xq*0 .+tol;
dam_x = [4701.18,4656.5]; 
dam_y = [4143.41,4392.1];
k = (dam_y[2]-dam_y[1])/(dam_x[2]-dam_x[1]);
flag = yq .- dam_y[1] - k*( xq .- dam_x[1] );
water_idx = findall(x->x<=0, flag);
h[water_idx] = 100 .- btm_q[water_idx] .+ tol;

bnd_y_idx = findall(y->abs( y .- 5250 ) < 150, yq);
bnd_x_idx = findall(x->abs( x .- 4500 ) < 200, xq);
bnd_idx = bnd_x_idx[findall(x->x in bnd_y_idx, bnd_x_idx)]
h[bnd_idx] .= tol;

h[:,ocean_idx] = tol .- btm_q[:,ocean_idx];

# btm_q = btm_q*0;
#test wb
# h = 110.0 .-btm_q;
# @show minimum(h)
if plot_init
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
    plot_data = Vp*Pq*(h+btm_q)

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
    fig1

    plt2 = Makie.mesh!(ax2, build_plot_data(plot_data, (rp,sp), (Vp*x, Vp*y)),
        color=vec(plot_data),
        # shading = false,
        colormap = :blues
        # colormap = :viridis
        );
    Makie.hidespines!(ax2)
    Makie.hidedecorations!(ax2)
    fig2
end
h0 = copy(h)
hu = h*0.0;
hv = h*0.0;

# t = 700
# folder = ""
# tag = string(folder, "malpasset_h");
# filename = string(tag, string(N),"_", string(Integer(round(t))), ".csv");
# h  = readdlm(filename, ',', Float64);
# tag = string(folder, "malpasset_hu");
# filename = string(tag, string(N),"_", string(Integer(round(t))), ".csv");
# hu  = readdlm(filename, ',', Float64);
# tag = string(folder, "malpasset_hv");
# filename = string(tag, string(N),"_", string(Integer(round(t))), ".csv");
# hv  = readdlm(filename, ',', Float64);
# error("end here")
"pack arguments into tuples"
Fmask = [findfirst(@. abs(rf[i] - rq) + abs(sf[i] - sq) < 100*eps()) for i in eachindex(rf)]
Nq = size(h,1); Nfq = size(E,1);nelem = size(h,2);
ops = ( Qrskew_ID,Qsskew_ID, Qr_ID, Qs_ID,
        Qrskew_ES,Qsskew_ES, QNr_sbp, QNs_sbp,
        E, M, M_inv, Mf_inv, Pf, Fmask);
dis_cst = (Cf, C, C_x, C_y, Nq, Nfq, nelem)
vgeo = (rxJ,sxJ,ryJ,syJ,J)
fgeo = (nxJ,nyJ,sJ, nx, ny)
nodemaps = (mapP,mapB)
(gQNxb, gQNyb) =  ESDG_bottom(QNr_sbp, QNs_sbp, btm_q, vgeo, g);
# U = (h, hu, hv, btm, gQNxb, gQNyb)
U = (h, hu, hv, btm_q)

rhs1_ID = zeros(Float64, size(h)); Mrhs1_ID = zeros(Float64, size(h));
rhs1_ES = zeros(Float64, size(h)); Mrhs1_ES = zeros(Float64, size(h));
lf = zeros(Float64, size(E*h));lff = zeros(Float64, size(h));
f1_ES = zeros(Float64, size(h)); f2_ES = zeros(Float64, size(h)); f3_ES= zeros(Float64, size(h));
f1_ID = zeros(Float64, size(h)); f2_ID = zeros(Float64, size(h)); f3_ID= zeros(Float64, size(h));
rhs1_CL = zeros(Float64, size(h)); rhs2_CL = zeros(Float64, size(h)); rhs3_CL = zeros(Float64, size(h));
u = zeros(Float64, size(h)); v = zeros(Float64, size(h));uf = zeros(Float64, size(E*h)); vf = zeros(Float64, size(E*h));
hf = zeros(Float64, size(E*h)); huf = zeros(Float64, size(E*h));hvf = zeros(Float64, size(E*h));
hP = zeros(Float64, size(E*h)); huP = zeros(Float64, size(E*h));hvP = zeros(Float64, size(E*h));
dh = zeros(Float64, size(E*h)); dhu = zeros(Float64, size(E*h));dhv = zeros(Float64, size(E*h));
lambdaf = zeros(Float64, size(E*h)); lambdaP = zeros(Float64, size(E*h)); cf = zeros(Float64, size(E*h));
h_L_next = zeros(Float64, size(h)); h_L_next_f = zeros(Float64, size(E*h))
h_H_next = zeros(Float64, size(h)); h_H_next_f = zeros(Float64, size(E*h))
f1_IDf = zeros(Float64, size(E*h)); f1_ESf = zeros(Float64, size(E*h));
rhs_1 = (zeros(Float64, size(h)), zeros(Float64, size(h)), zeros(Float64, size(h)))
rhs_2 = (zeros(Float64, size(h)), zeros(Float64, size(h)), zeros(Float64, size(h)))
htmp = zeros(Float64, size(h)); hutmp = zeros(Float64, size(h)); hvtmp = zeros(Float64, size(h));
FS1 = zeros(Float64, size(E*h)); FS2 = zeros(Float64, size(E*h)); FS3 = zeros(Float64, size(E*h));
dt_local =  zeros(Float64, (1, size(h,2))); L_E =  zeros(Float64, (1, size(h,2)))
f1_ES_e = zeros(Float64, (Nq,Nq));f2_ES_e = zeros(Float64, (Nq,Nq));f3_ES_e = zeros(Float64, (Nq,Nq))
f1_ID_e = zeros(Float64, (Nq,Nq));f2_ID_e = zeros(Float64, (Nq,Nq));f3_ID_e = zeros(Float64, (Nq,Nq))
UL_E = (Float64, Float64, Float64); UR_E = (Float64, Float64, Float64); vgeo_e =(Float64, Float64, Float64, Float64);
pre_allo = (u, v, uf, vf, hf, huf, hvf, hP, huP, hvP, dh, dhu, dhv,
            rhs1_ID, Mrhs1_ID, rhs1_ES, Mrhs1_ES, lf, lff, f1_ES, f2_ES, f3_ES, f1_ID, f2_ID, f3_ID,
            rhs1_CL, rhs2_CL, rhs3_CL, lambdaf, lambdaP, cf, UL_E, UR_E, vgeo_e,
            h_L_next, h_L_next_f, f1_IDf, f1_ESf, FS1, FS2, FS3, dt_local, L_E)