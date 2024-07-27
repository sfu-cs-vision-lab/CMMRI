
% for details on the terminology see the following publications:
% "An object-color space", A.D. Logvinenko, Journal of Vision
% "Metamer Mismatch Volumes", Logvinenko, Godau and Funt, Proc. CGIV 2012


classdef ColorMechanism < handle
    properties 
        illuminant
        cmf
        wavelengths
        cum_cmf_pp
        white
        n_sensors
        invomegaspline
        omegaspline
    end
    
    methods
        % constructor
        % each color mechanism consists of illuminant and cmf, which need
        % to be provided on creation, along with the corresponding
        % wavelengths
        % illuminant normalized to make Y (second column of cmf) = 100
        function obj = ColorMechanism(illuminant, cmf, wavelengths,no_norm_ill)
            obj.illuminant=illuminant;
            obj.cmf=cmf;
            obj.wavelengths=wavelengths;
            if nargin<4
                obj.normalize_illuminant();
            end
            obj.n_sensors=size(cmf,2);
            obj.generate_interpolator();
        end
        
        function [ ill_norm ] = normalize_illuminant(obj, normvalue, normchannel)
            
            % normalizes the cmf to that maximum of each channel is 100
            n_channels=size(obj.cmf,2);
            if nargin<2
                normvalue=100;
            end
            if nargin<3
                normchannel=ceil(n_channels/2);
            end
            
            max_rgb=obj.illuminant(:)'*repmat(obj.cmf(:,normchannel),[1 n_channels]);
            obj.cmf=obj.cmf./repmat(max_rgb,size(obj.cmf,1),1)*normvalue;
            obj.white=obj.illuminant'*obj.cmf;
        end
        
        % generates the inperpolation functions needed to quickly calculate
        % responses        
        function generate_interpolator(obj)
            % create interpolation functions for sensor responses
            illcones=obj.cmf.*repmat(obj.illuminant,1,obj.n_sensors);
            cumcones=cumtrapz(illcones); % cumulative cones
            obj.white=cumcones(end,:); % maximal response (white)
            obj.white=obj.white(:)';
            cumcones=repmat(cumcones(end,:),size(cumcones,1),1)-cumcones; % invert the response (integrate from wl to end instead of integrating from start to wl)
            obj.cum_cmf_pp=spline(obj.wavelengths',cumcones'); % interpolation polygon
            
            % create the de/reeparametrization function interpolators
            sigma=cumtrapz(sqrt(sum(obj.cmf.^2,2)));
            sigma=sigma/max(sigma);
            [~,idx]=unique(sigma(end:-1:1,:));
            idx=length(sigma)-idx+1;
            obj.invomegaspline=spline(sigma(idx), obj.wavelengths(idx));
            obj.omegaspline=spline(obj.wavelengths(idx), sigma(idx));
        end
        % deparamatrization: from [0,1] to wavelengths
        function out=deparametrize(obj, wl)
            out=ppval(obj.invomegaspline,wl);
        end
        
        % reparametrization: from wavelengths to [0,1]
        function out=reparametrize(obj, wl)
            out=ppval(obj.omegaspline,wl);
        end
        
        % this will evaluate the interpolation polygon at the given
        % wavelengths
        % corresponds to integration from lambda to lambda_max
        function out=getcum(obj, wavelengths)
            out=ppval(obj.cum_cmf_pp, wavelengths);
        end
        
        % convert an n-transition optimal spectrum to sensor responses
        % each line in ntrs is a spectrum with l1<l2<l3... for type 1 or
        % l1>l2>l3... for type 2
        function out=ntr2lms(obj,ntrs, types)
            % for more then one transition wavelength, the order determines
            % the type, otherwise the types argument needs to be given
            % set to 1 for type 2 reflectances
            if ~exist('types','var')
                types=(ntrs(:,1)>ntrs(:,end))+1;
            end
            out=zeros(size(ntrs,1),obj.n_sensors);
            for i=1:size(ntrs,1)
                current_ntr=ntrs(i,:);
                v=obj.getcum(sort(current_ntr));
                lms=sum([v(:,1:2:end) -v(:,2:2:end)],2)';
                if types(i)==2
                    % type2 response is white-type1
                    lms=obj.white-lms;
                end
                out(i,:)=lms;
            end
        end
        % wrapper function for reparametrized optimal wl to lms
        % conversion
        % a _r suffic indicates reparametrization
        function out=ntr2lms_r(obj,ntrs_r, types)
%             if nargin < 2
%                 % set to 1 for type 2 reflectances
%                 types=ntrs_r(:,1)>ntrs_r(:,end);
%             end
            if ~exist('types')
                types=ntrs(:,1)>ntrs(:,end)+1;
            end
            ntrs=obj.deparametrize(ntrs_r);
            out=obj.ntr2lms(ntrs, types);
        end
        
        function [ LMS ] = rect2lms( obj, rect_wl )
            % convert a rectangular metamer given by transition wavelengths
            % and purity to lms
            % rect_wl contains columns for purity, l_1, l_2
            
            % calculate LMS for optimal reflectance function given by
            % transition wavelengths
            rect_lms=obj.ntr2lms(rect_wl(:,2:3));

            % calculate final LMS using given alpha as a convex combination
            % of oLMS and grey
            nSamples=size(rect_wl,1);
            grey=repmat(obj.white/2,nSamples,1);
            opt_lms_nogrey=rect_lms-grey; % subtract grey from optimal lms
            LMS=grey+repmat(rect_wl(:,1),1,3).*opt_lms_nogrey; % add optimal lms weighted by purity
        end
        
        function [ rect ] = adl2opt( obj, adl )
            % converts the provided adl to a rectangular metamer
            % rect has 3 columns: purity + 2 transition wavelengths
            % purity is not changed by this function
            % type is determined by the order of the transition wavelengths

            rect=adl;
            rect(:,2:3)=0;

            wlrange=obj.wavelengths(end)-obj.wavelengths(1);
            l1=adl(:,3)-adl(:,2)/2;
            l2=adl(:,3)+adl(:,2)/2;
            % correct for type 2
            l1(l1<obj.wavelengths(1))=l1(l1<obj.wavelengths(1))+wlrange;
            l2(obj.wavelengths(end)<l2)=l2(obj.wavelengths(end)<l2)-wlrange;
            rect(:,2:3)=[l1 l2];
        end
        
        function out=adl2lms(obj, adl)
            rect=obj.adl2opt(adl);
            out=obj.rect2lms(rect);
        end
        
        function [vol XYZ]=draw_objcolorsolid(obj, alpha)
            % draws a colored 3d model of the object color solid in XYZ
            % space into the current figure
            % this obviously only works for 3 dimensions, and assumes the
            % values are something like XYZ to draw colors
            if nargin<2
                alpha=0.5;
            end

            % generate samples on ocs surface using reparametrized values for a more even
            % distribution of points

                range=linspace(0,1,100);
                N=numel(range)^2;

                rOPT=zeros(N,2);

                i=0;
                for L1=range
                    for L2=range
                        i=i+1;
                        rOPT(i,:)=[L1,L2];
                    end
                end

                % convert to normal wl in nm
                OPT=obj.deparametrize(rOPT);

            % convert to XYZ

            XYZ=obj.ntr2lms(OPT);

            % draw the convex hull

            vol=draw_conv_hull_color(XYZ,alpha)
        end
        

        
        function [ adl ] = opt2adl( obj, rect, wl )
            % converts the provided rectangular metamers to adl
            % rect has 3 columns: purity + 2 transition wavelengths
            % purity is not changed by this function
            % type is determined by the order of the transition wavelengths

            adl=rect;
            adl(:,2:3)=0;
            %calculate lambda central
            lc=mean(rect(:,2:3),2);
            % correction for type 2
            t2=rect(:,3)<rect(:,2);
            wlrange=obj.wavelengths(end)-obj.wavelengths(1);
            lc(t2)=lc(t2)-wlrange/2;
            lc(lc<wl(1))=lc(lc<wl(1))+wlrange;
            adl(:,3)=lc;

            %calculate delta
            d=abs(rect(:,2)-rect(:,3));
            d(t2)=wlrange-d(t2);
            adl(:,2)=d;
        end
        
        function [hull_x, hull_y, hull_l]=draw_ocs_2d(obj)
            % draw the objects color solid in 2d xy chromaticity space
            n_points=2000; % how many points on the boundary?
            wls=linspace(obj.wavelengths(1)+0.5, obj.wavelengths(end)-0.5, n_points);
            step_size=obj.wavelengths(2)-obj.wavelengths(1);
            % make a list of "pure" colors
            ntrs=[wls(:)-step_size/2 wls(:)+step_size/2];
            lms=obj.ntr2lms(ntrs);
            
            x= lms(:,1)./sum(lms,2);
            y= lms(:,2)./sum(lms,2);
            l= lms(:,2);
            
            k = convhull(x,y);
            plot(x(k),y(k),'-');
            hull_x=x(k);
            hull_y=y(k);
            hull_l=l(k);
        end
        
        function [vol]=get_ocs_vol(obj)
            % get the approximate volume of the object color solid
            
            
            % generate samples on ocs surface using reparametrized values for a more even
            % distribution of points
            
                range=linspace(0,1,100);
                N=numel(range)^2;

                rOPT=zeros(N,2);

                i=0;
                for L1=range
                    for L2=range
                        i=i+1;
                        rOPT(i,:)=[L1,L2];
                    end
                end

                % convert to normal wl in nm
                OPT=obj.deparametrize(rOPT);

            % convert to XYZ

            XYZ=obj.ntr2lms(OPT);

            % draw the convex hull

            [~,vol]=convhull(XYZ);

        end
        
    end
end