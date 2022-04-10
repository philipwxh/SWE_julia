global tau = 1
avg(a,b) = .5*(a+b)
@inline function fS2D(UL,UR,g)
    hL,huL,hvL = UL
    hR,huR,hvR = UR
    uL,vL = (x->x./hL).((huL,hvL))
    uR,vR = (x->x./hR).((huR,hvR))
    fxS1 = @. avg(huL,huR)
    fxS2 = @. avg(huL,huR)*avg(uL,uR) + .5*g*hL*hR
    fxS3 = @. avg(huL,huR)*avg(vL,vR)

    fyS1 = @. avg(hvL,hvR)
    fyS2 = @. avg(hvL,hvR)*avg(uL,uR)
    fyS3 = @. avg(hvL,hvR)*avg(vL,vR) + .5*g*hL*hR
    return (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3)
end

function fS2D_LF(UL,UR,g)
    hL,huL,hvL = UL
    hR,huR,hvR = UR
    uL,vL = (x->x./hL).((huL,hvL))
    uR,vR = (x->x./hR).((huR,hvR))
    fxS1 = @. avg(huL,huR)
    fxS2 = @. avg(huL*uL,huR*uR) + .5*g*avg(hL*hL,hR*hR)
    fxS3 = @. avg(huL*vL,huR*vR)

    fyS1 = @. avg(hvL,hvR)
    fyS2 = @. avg(hvL*uL,hvR*uR)
    fyS3 = @. avg(hvL*vL,hvR*vR) + .5*g*avg(hL*hL,hR*hR)
    return (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3)
end
function swe_2d_esdg_surface(UL, UR, dU, Pf, fgeo, c, g)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}
    nxJ,nyJ,sJ,nx,ny = fgeo
    (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D(UL,UR,g)
    dh, dhu, dhv = dU
    (hf, huf, hvf)=UL;
    (hP, huP, hvP)=UR;
    fs1 = @. fxS1*nxJ + fyS1*nyJ;
    fs2 = @. fxS2*nxJ + fyS2*nyJ;
    fs3 = @. fxS3*nxJ + fyS3*nyJ;
    f1_ES =  Pf*(fs1 .- .5*tau*c.*(hP.-hf).*sJ);
    f2_ES =  Pf*(fs2 .- .5*tau*c.*(huP.-huf).*sJ);
    f3_ES =  Pf*(fs3 .- .5*tau*c.*(hvP.-hvf).*sJ);
    return f1_ES, f2_ES, f3_ES
end

function swe_2d_esdg_n_surface(UL, UR, dU, fgeo, g)::Tuple{Float64,Float64,Float64}
    nxJ,nyJ,sJ,c = fgeo
    (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D_LF(UL,UR,g)
    # (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D(UL,UR,g)
    dh, dhu, dhv = dU
    fs1 = fxS1*nxJ + fyS1*nyJ;
    fs2 = fxS2*nxJ + fyS2*nyJ;
    fs3 = fxS3*nxJ + fyS3*nyJ;
    fs1 = fs1 - .5*tau*c*dh*sJ;
    fs2 = fs2 - .5*tau*c*dhu*sJ
    fs3 = fs3 - .5*tau*c*dhv*sJ
    # f1_ID = Pf * (0.5*fs1 .- .5*tau*c.*dh.*sJ)
    # f2_ID = Pf * (0.5*fs2 .- .5*tau*c.*dhu.*sJ)
    # f3_ID = Pf * (0.5*fs3 .- .5*tau*c.*dhv.*sJ)
    return fs1, fs2, fs3
end

function swe_2d_esdg_vol(UL_E, UR_E, ops, vgeo_e, i, j, btm, g)
    Qr_ID,Qs_ID,Qrb_ID,Qsb_ID,Qr_ES,Qs_ES,QNr_sbp, QNs_sbp, E,M_inv,Pf= ops
    (rxJ_i, sxJ_i, ryJ_i, syJ_i) = vgeo_e;
    (FxV1,FxV2,FxV3),(FyV1,FyV2,FyV3) = fS2D(UL_E,UR_E,g)
    h_i, hu_i, hv_i = UL_E; h_j, hu_j, hv_j = UR_E;
    # QNx_ij = Qr_ES[i,j]*(rxJ[i,e]+ rxJ[j,e]) + Qs_ES[i,j]*(sxJ[i,e]+ sxJ[j,e]);
    # QNy_ij = Qr_ES[i,j]*(ryJ[i,e]+ ryJ[j,e]) + Qs_ES[i,j]*(syJ[i,e]+ syJ[j,e]);
    QNx_ij = Qr_ES[i,j]*(rxJ_i*2) + Qs_ES[i,j]*(sxJ_i*2);
    QNy_ij = Qr_ES[i,j]*(ryJ_i*2) + Qs_ES[i,j]*(syJ_i*2);

    QNxb_ij = QNr_sbp[i,j]*(rxJ_i) + QNs_sbp[i,j]*(sxJ_i);
    QNyb_ij = QNr_sbp[i,j]*(ryJ_i) + QNs_sbp[i,j]*(syJ_i);

    fv1_ES = (QNx_ij*FxV1 + QNy_ij*FyV1);
    # fv2_ES = (QNx_ij*FxV2 + QNy_ij*FyV2);
    # fv3_ES = (QNx_ij*FxV3 + QNy_ij*FyV3);

    # fv2_i_ES = (QNx_ij*FxV2 + QNy_ij*FyV2);
    # fv3_i_ES = (QNx_ij*FxV3 + QNy_ij*FyV3);
    # fv2_j_ES = -(QNx_ij*FxV2 + QNy_ij*FyV2);
    # fv3_j_ES = -(QNx_ij*FxV3 + QNy_ij*FyV3);
    fv2_i_ES = (QNx_ij*FxV2 + QNy_ij*FyV2) + g*h_i*QNxb_ij*btm[j];
    fv3_i_ES = (QNx_ij*FxV3 + QNy_ij*FyV3) + g*h_i*QNyb_ij*btm[j];
    fv2_j_ES = -(QNx_ij*FxV2 + QNy_ij*FyV2) - g*h_j*QNxb_ij*btm[i];
    fv3_j_ES = -(QNx_ij*FxV3 + QNy_ij*FyV3) - g*h_j*QNyb_ij*btm[i];
    return fv1_ES, fv2_i_ES, fv3_i_ES, -fv1_ES, fv2_j_ES, fv3_j_ES
end

function swe_2d_ES_h(UL_E,UR_E, Qr_ES, Qs_ES, vgeo_e, i, j, g)
    (rxJ_i, sxJ_i, ryJ_i, syJ_i) = vgeo_e;
    (FxV1,FxV2,FxV3),(FyV1,FyV2,FyV3) = fS2D(UL_E,UR_E,g)
    # h_i, hu_i, hv_i = UL_E; h_j, hu_j, hv_j = UR_E;
    # QNx_ij = Qr_ES[i,j]*(rxJ[i,e]+ rxJ[j,e]) + Qs_ES[i,j]*(sxJ[i,e]+ sxJ[j,e]);
    # QNy_ij = Qr_ES[i,j]*(ryJ[i,e]+ ryJ[j,e]) + Qs_ES[i,j]*(syJ[i,e]+ syJ[j,e]);
    QNx_ij = Qr_ES[i,j]*(rxJ_i*2) + Qs_ES[i,j]*(sxJ_i*2);
    QNy_ij = Qr_ES[i,j]*(ryJ_i*2) + Qs_ES[i,j]*(syJ_i*2);


    fv1_ES = (QNx_ij*FxV1 + QNy_ij*FyV1);
    return fv1_ES, -fv1_ES
end

function swe_2d_ID_surface(UL, UR, dU, Pf, fgeo, c, g)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}
    nxJ,nyJ,sJ,nx,ny = fgeo
    (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D_LF(UL,UR,g)
    # (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D(UL,UR,g)
    dh, dhu, dhv = dU
    fs1 = @. fxS1*nxJ + fyS1*nyJ;
    fs2 = @. fxS2*nxJ + fyS2*nyJ;
    fs3 = @. fxS3*nxJ + fyS3*nyJ;
    f1_ID = Pf * (0.5*fs1 .- .5*tau*c.*dh.*sJ)
    f2_ID = Pf * (0.5*fs2 .- .5*tau*c.*dhu.*sJ)
    f3_ID = Pf * (0.5*fs3 .- .5*tau*c.*dhv.*sJ)
    # f1_ID = Pf * (0.5*fs1) - transpose(E)*(Cf.*c.*dh)
    # f2_ID = Pf * (0.5*fs2) - transpose(E)*(Cf.*c.*dhu)
    # f3_ID = Pf * (0.5*fs3) - transpose(E)*(Cf.*c.*dhv)
    return f1_ID, f2_ID, f3_ID
end

function swe_2d_ID_n_surface(UL, UR, dU, fgeo, g)::Tuple{Float64,Float64,Float64}
    nxJ,nyJ,sJ,c = fgeo
    (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D_LF(UL,UR,g)
    # (fxS1,fxS2,fxS3),(fyS1,fyS2,fyS3) = fS2D(UL,UR,g)
    dh, dhu, dhv = dU
    fs1 = fxS1*nxJ + fyS1*nyJ;
    fs2 = fxS2*nxJ + fyS2*nyJ;
    fs3 = fxS3*nxJ + fyS3*nyJ;
    fs1 = 0.5*fs1 - .5*tau*c*dh*sJ;
    fs2 = 0.5*fs2 - .5*tau*c*dhu*sJ
    fs3 = 0.5*fs3 - .5*tau*c*dhv*sJ
    # f1_ID = Pf * (0.5*fs1) - transpose(E)*(Cf.*c.*dh)
    # f2_ID = Pf * (0.5*fs2) - transpose(E)*(Cf.*c.*dhu)
    # f3_ID = Pf * (0.5*fs3) - transpose(E)*(Cf.*c.*dhv)
    return fs1, fs2, fs3
end

function swe_2d_ID_vol(UL_E, ops, vgeo_e, i, j, btm, g)
    Qr_ID,Qs_ID,Qrb_ID,Qsb_ID,Qr_ES,Qs_ES,QNr_sbp, QNs_sbp, E,M_inv,Pf= ops
    (rxJ_i, sxJ_i, ryJ_i, syJ_i) = vgeo_e;
    (fxV1,fxV2,fxV3),(fyV1,fyV2,fyV3) = fS2D_LF(UL_E,UL_E,g)
    # (fxV1,fxV2,fxV3),(fyV1,fyV2,fyV3) = fS2D(UL_E,UL_E,g)
    h_i, hu_i, hv_i = UL_E;

    QNx_ij = 2*(Qr_ID[i,j]*rxJ_i + Qs_ID[i,j]*sxJ_i);
    QNy_ij = 2*(Qr_ID[i,j]*ryJ_i + Qs_ID[i,j]*syJ_i);

    QNxb_ij = Qrb_ID[i,j]*rxJ_i + Qsb_ID[i,j]*sxJ_i;
    QNyb_ij = Qrb_ID[i,j]*ryJ_i + Qsb_ID[i,j]*syJ_i;

    # QNxb_ij = QNr_sbp[i,j]*(rxJ_i) + QNs_sbp[i,j]*(sxJ_i);
    # QNyb_ij = QNr_sbp[i,j]*(ryJ_i) + QNs_sbp[i,j]*(syJ_i);

    fv1_ID = QNx_ij*fxV1 + QNy_ij*fyV1;
    fv2_ID = QNx_ij*fxV2 + QNy_ij*fyV2 + g*h_i*QNxb_ij*btm[j];
    fv3_ID = QNx_ij*fxV3 + QNy_ij*fyV3 + g*h_i*QNyb_ij*btm[j];
    return fv1_ID, fv2_ID, fv3_ID
end

function swe_2d_ID_h(UL_E, Qr_ID, Qs_ID, vgeo_e, i, j, g)
    (rxJ_i, sxJ_i, ryJ_i, syJ_i) = vgeo_e;
    (fxV1,fxV2,fxV3),(fyV1,fyV2,fyV3) = fS2D_LF(UL_E,UL_E,g)
    # (fxV1,fxV2,fxV3),(fyV1,fyV2,fyV3) = fS2D(UL_E,UL_E,g)
    Qr_ID_ij = Qr_ID[i,j]; Qs_ID_ij = Qs_ID[i,j];
    dhdx  = rxJ_i*(Qr_ID_ij*fxV1) + sxJ_i*(Qs_ID_ij*fxV1)
    dhdy  = ryJ_i*(Qr_ID_ij*fyV1) + syJ_i*(Qs_ID_ij*fyV1)
    fv1_ID = 2*(dhdx  + dhdy)
    return fv1_ID
end

# function zero_vel!(h, hu, hv, tol, val = 0)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}
#     (m,n) = size(h);
#     for i = 1:m
#         for j = 1:n
#             if h[i,j]< 100*tol
#                 # h[i,j] = tol;
#                 hu[i,j] = val;
#                 hv[i,j] = val;
#             end
#             # if abs(hu[i,j]/h[i,j])>50
#             #     hu[i,j] = sign(hu[i,j]) * 50*h[i,j]
#             # end
#             # if abs(hv[i,j]/h[i,j])>50
#             #     hv[i,j] = sign(hv[i,j]) * 50*h[i,j]
#             # end
#             # if h[i,j]>1.1*50
#             #     h[i,j] = 50*1.1
#             # end
#             if h[i,j]<tol/2
#                 h[i,j] = tol;
#             end
#         end
#     end
#     return h, hu, hv
# end

function zero_vel!(h, hu, hv, tol, h0, ocean_idx, val = 0)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}
    (m,n) = size(h);
    for e = 1:n
        # if minimum(h[:,e]) < 0
        #     h[:,e] .= mean(h[:,e]);
        #     hu[:,e] .=0
        #     hv[:,e] .=0
        # end
    	# k = rand()
        for i = 1:m
            if h[i,e] < tol #&& h[i,e] > -tol*10
                h[i,e] = tol
            end
            if h[i,e]< 10*tol
                # h[i,e] = tol;
                hu[i,e] = val;
                hv[i,e] = val;
            end
            # if abs(hu[i,j]/h[i,j])> 20 + k*2
            #     hu[i,j] = sign(hu[i,j]) * (20 + k*2)*h[i,j]
            # end
            # if abs(hv[i,j]/h[i,j])> 10 + k
            #     hv[i,j] = sign(hv[i,j]) * (10 + k)*h[i,j]
            # end
            # if h[i,j]>h_max - (10 + k*8)*tol
            #     h[i,j] = h_max - (10 + k*8)*tol
            # end
            # if h[i,j]<tol/2
            #     h[i,j] = tol;
            # end
        end
    end
    h[ocean_idx] = h0[ocean_idx]
    hu[ocean_idx] .= 0
    hv[ocean_idx] .= 0
    return h, hu, hv
end

function vel_reg!(h, hu, hv, h_0_max, tol, val = 0)::Tuple{Array{Float64,2},Array{Float64,2},Array{Float64,2}}
    (m,n) = size(h);
    for i = 1:m
        for j = 1:n
            h_e = max(h[i,j], tol, h_0_max);
            if h[i,j]<tol*h_0_max
                h[i,j] = tol;
                hu[i,j] = 2*h[i,j]/(h[i,j]^2+h_e^2)*hu[i,j];
                hv[i,j] = 2*h[i,j]/(h[i,j]^2+h_e^2)*hv[i,j];
            end
            # if h[i,j]<tol/2
            #     # h[i,j] = tol
            #     h[:,j] .= tol
            # end
        end
    end
    return h, hu, hv
end

function set_val_com!(A, B, tol, val = 0)::Array{Float64,2}
    (m,n) = size(A);
    for i = 1:m
        for j = 1:n
            if A[i,j]<tol
                B[i,j] = val;
            end
        end
    end
    return B
end

function SWE_vortex( x, y, t )
    # This function calculate the analytical solution of vortex translate
    # across the domain without changing its shape .
    # H = H_inf - ?^2 / (32 * pi^2) * e^( -2(r^2 - 1) ),
    # u = u_inf - ?   / (2 * pi)    * e^(  -(r^2 - 1) ) * yt
    # v = v_inf + ?   / (2 * pi)    * e^(  -(r^2 - 1) ) * xt,
    H_inf = 1;
    u_inf = 1;
    v_inf = 0;
    xc = 0;
    yc = 0;
    xt = x .- xc .- u_inf*t;
    yt = y .- yc .- v_inf*t;
    r_sq = xt.^2 + yt.^2;
    beta = 5;
    g = 2;
    H = H_inf .- beta^2 / (32 * pi^2) * exp.( -2 * (r_sq .- 1) );
    u = u_inf .- beta   / (2 * pi)    * exp.(  -( r_sq .- 1) ) .* yt;
    v = v_inf .+ beta   / (2 * pi)    * exp.(  -( r_sq .- 1) ) .* xt;
    return H, u, v
end