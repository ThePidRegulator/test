
packs = pkg('list');
for jj = 1:numel(packs),
        pkg('load', packs{jj}.name);
end