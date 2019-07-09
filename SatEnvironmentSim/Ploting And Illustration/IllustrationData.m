classdef IllustrationData < handle
        properties (SetAcces = public, GetAcces = public)
                S = [];
                rSun = {};
                rSat = {};
                albedoDir = {};
                albedoIrr = {};
                albedoSrce = {};
                albedoMaxIrr = [];
                satFovMap = {};
                sunlitMap = {};
                litFovMap = {};
                reflMap = [];
                EarthH = [];
        end
end
