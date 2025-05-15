function d_3d_matrix = fd_for_3d_matrix(M,dx,dy,type)

d_3d_matrix = zeros(size(M));
switch type
    case 'ddx_fwd'
        d_3d_matrix(1,:,:) = ddx_fwd(squeeze(M(1,:,:)),dx);
        d_3d_matrix(2,:,:) = ddx_fwd(squeeze(M(2,:,:)),dx);
        d_3d_matrix(3,:,:) = ddx_fwd(squeeze(M(3,:,:)),dx);
        d_3d_matrix(4,:,:) = ddx_fwd(squeeze(M(4,:,:)),dx);
    case 'ddx_bwd'
        d_3d_matrix(1,:,:) = ddx_bwd(squeeze(M(1,:,:)),dx);
        d_3d_matrix(2,:,:) = ddx_bwd(squeeze(M(2,:,:)),dx);
        d_3d_matrix(3,:,:) = ddx_bwd(squeeze(M(3,:,:)),dx);
        d_3d_matrix(4,:,:) = ddx_bwd(squeeze(M(4,:,:)),dx);
    case 'ddy_bwd'
        d_3d_matrix(1,:,:) = ddy_bwd(squeeze(M(1,:,:)),dy);
        d_3d_matrix(2,:,:) = ddy_bwd(squeeze(M(2,:,:)),dy);
        d_3d_matrix(3,:,:) = ddy_bwd(squeeze(M(3,:,:)),dy);
        d_3d_matrix(4,:,:) = ddy_bwd(squeeze(M(4,:,:)),dy);
    case 'ddy_fwd'
        d_3d_matrix(1,:,:) = ddy_fwd(squeeze(M(1,:,:)),dy);
        d_3d_matrix(2,:,:) = ddy_fwd(squeeze(M(2,:,:)),dy);
        d_3d_matrix(3,:,:) = ddy_fwd(squeeze(M(3,:,:)),dy);
        d_3d_matrix(4,:,:) = ddy_fwd(squeeze(M(4,:,:)),dy);

    otherwise
        disp('Warning')

end
end
