function SliceInfo_out = yGetSliceInfo(SliceInfo)

Nslice = length(SliceInfo);
SliceInfo_out.SliceThickness = zeros(Nslice,1,'single');
SliceInfo_out.FOV = zeros(Nslice,2,'single');
SliceInfo_out.InPlaneRot = zeros(Nslice,1,'single');
SliceInfo_out.SlicePosition = zeros(Nslice,3,'single');
SliceInfo_out.SliceOrintation = zeros(Nslice,3,'single');

for Nsl=1:length(SliceInfo)
    SliceInfo_out.SliceThickness(Nsl) = SliceInfo(Nsl).dThickness;
    SliceInfo_out.FOV(Nsl,1) = SliceInfo(Nsl).dPhaseFOV;
    SliceInfo_out.FOV(Nsl,2) = SliceInfo(Nsl).dReadoutFOV;
    if isfield(SliceInfo(Nsl),'dInPlaneRot')
        if ~isempty(SliceInfo(Nsl).dInPlaneRot)
            SliceInfo_out.InPlaneRot(Nsl) = SliceInfo(Nsl).dInPlaneRot;
        end
    end
    if isfield(SliceInfo(Nsl).sPosition,'dSag')
        SliceInfo_out.SlicePosition(Nsl,1) = SliceInfo(Nsl).sPosition.dSag;
    end
    if isfield(SliceInfo(Nsl).sPosition,'dCor')
        SliceInfo_out.SlicePosition(Nsl,2) = SliceInfo(Nsl).sPosition.dCor;
    end
    if isfield(SliceInfo(Nsl).sPosition,'dTra')
        SliceInfo_out.SlicePosition(Nsl,3) = SliceInfo(Nsl).sPosition.dTra;
    end
    if isfield(SliceInfo(Nsl).sNormal,'dSag')
        SliceInfo_out.SliceOrintation(Nsl,1) = SliceInfo(Nsl).sNormal.dSag;
    end
    if isfield(SliceInfo(Nsl).sNormal,'dCor')
        SliceInfo_out.SliceOrintation(Nsl,2) = SliceInfo(Nsl).sNormal.dCor;
    end
    if isfield(SliceInfo(Nsl).sNormal,'dTra')
        SliceInfo_out.SliceOrintation(Nsl,3) = SliceInfo(Nsl).sNormal.dTra;
    end
end

end