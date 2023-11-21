classdef loaders_class<handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dirpath
        hullpath
        movpath
        segpath
        hullfile
        hullRecfile
        segfile
        easypath
        mov
        hull3dname
        struct_hullfile
    end
    
    methods
        function obj = loaders_class(path,mov,easyname,varargin)
            parser = inputParser;
            addParameter(parser,'hullfile','\\hull_op\'); % number of camera pointing in Z lab axis
            addParameter(parser,'segfile','\\Segmentation\'); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            hullfile = parser.Results.hullfile;
            segfile = parser.Results.segfile;
           
            % initilize loaders class
            obj.mov = mov;
            movdir = sprintf('\\mov%d',mov);
            obj.hullfile  = sprintf('hull_mov%d',mov);
            obj.struct_hullfile  = sprintf('Shull_mov%d',mov);
            obj.hullRecfile  = sprintf('hullRec_mov%d',mov);
            obj.segfile  = sprintf('mov%d_seg',mov);
            obj.dirpath = path;
            obj.movpath = [path,movdir];
            obj.easypath  = [path,'\',easyname];
            obj.hullpath  = [path,movdir,hullfile];
            obj.segpath  = [path,movdir,segfile];
            obj.hull3dname = sprintf('hull3d_mov%d',mov);
        end
        
        function easy = easywand(obj)
            % load easywand file
            a = load(obj.easypath);
            easy = a.easyWandData;
        end
        
        function sp = loadsparse(obj)
            % load movie sparse files
            listing = dir([obj.movpath]);
            spcount = 1;
            for k = 1:1:length(listing)
                if strfind(listing(k).name,'sparse')
                    sp{spcount} = load([obj.movpath,'\',listing(k).name]);
                    spcount = spcount + 1;
                end
            end
        end
        
        function hull_op = loadhullRec(obj,loadflg)
            % load hull_op
            hull_op = [];
            if exist(obj.hullpath,'dir') == 0
                mkdir(obj.hullpath);
            end
            if loadflg == 1
                filename = sprintf('hullRec_mov%d',obj.mov);
                a = load([obj.hullpath,'\',filename]);
                hull_op = a.hullRec;
            end
        end
        
        function [hull,hull3d] = loadhull(obj,loadflg,varargin)
            % load hull
            parser = inputParser;
            addParameter(parser,'load_hull3d',1); % number of camera pointing in Z lab axis
            parse(parser, varargin{:})
            hull = [];hull3d = [];
            if loadflg == 1
                filename = sprintf('hull_mov%d',obj.mov);
                a = load([obj.hullpath,'\',filename]);
                hull = a.hull;
                if parser.Results.load_hull3d == 1
                    b = load([obj.hullpath,obj.hull3dname],'hull3d');
                    hull3d = b.hull3d;
                end
            end
        end
        
        function seg = loadSegfile(obj,loadflg)
            seg = [];
            if exist(obj.segpath,'dir') == 0
                mkdir(obj.segpath);
            end
            if loadflg == 1
                a = load([obj.segpath,'\',obj.segfile]);
                seg = a.seg;
            end
        end
    end
end

