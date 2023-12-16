function out = cx_read_tif_tzyx(filename_tif, ...
    out_datatype, size_read_to, ...
    size_z_read_from, size_t_read_from, ...
    inds_z_read_from, inds_t_read_from)

allow_data_type_conversion = 0;

%scanimage writes as tzyx, and so does my python registration and roi extraction
%so read it like this, then permute output later (faster that way)

tr = Tiff(filename_tif, 'r');

countz_ti = 0;
for ti =  1:size_t_read_from
    if  ismember(ti, inds_t_read_from)
        countz_ti = countz_ti + 1;
    end
    countz_zi = 0;
    % if mod(ti, 100)==0
    %     display(["on frame " num2str(ti)])
    % end
    for zi = 1:size_z_read_from
        if ismember(zi, inds_z_read_from) & ismember(ti, inds_t_read_from)
            countz_zi = countz_zi + 1;
            if countz_ti==1 & countz_zi==1 %on first frame, find data type and create output array
                tmpframe = tr.read();
                in_datatype = class(tmpframe);
                out = zeros(size_read_to, in_datatype);
                out(countz_ti,countz_zi,:,:) = tmpframe; 
            else
                out(countz_ti,countz_zi,:,:) = tr.read(); %tiff read faster than imread
            end
        end
        try
            tr.nextDirectory()
        catch
            "FINAL TIF FRAME"
        end
    end
end
%
% out = permute(out, [2 3 1]);
% out = reshape(out, size_read_to(2), size_read_to(3), length(inds_t_read_from), length(inds_z_read_from));
out = permute(out, [3 4 2 1]); %reshape into y x z t

%make minimal required adjustments to switch data types
if ~strcmp(in_datatype, out_datatype)
    if ~allow_data_type_conversion
        "ERROR, DATA TYPE MISMATCH"
        error
    else
        "WARNING, DATA TYPE MISMATCH, CONVERTING"
        switch in_datatype
            case'double'
                error

            case 'single'
                error

            otherwise
                inprec = str2double(regexp(in_datatype,'\d*','Match'));
                outprec = str2double(regexp(out_datatype,'\d*','Match'));
                minout = min(out(:));
                if minout<0 & startsWith(out_datatype, 'uint')
                    out = single(out);
                    out = out - double(minout);
                end
                if inprec<=outprec
                    eval(['out = ' out_datatype '(out);'])
                elseif inprec>outprec
                    maxout = max(out(:));
                    if maxout > 2^outprec-1
                        "ERROR, CLIPPING REQUIRED, CHANGE OUTPUT TYPE"
                        error
                    else
                        eval(['out = ' out_datatype '(out);'])
                    end
                end
        end
    end
end

close(tr)

