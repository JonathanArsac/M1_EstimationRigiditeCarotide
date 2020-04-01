function map_hist = get_map_hist(I)
sum_hist = 0;
hist_cumul = zeros(1, 255);
for i = 1:255
    sum_hist = sum_hist + sum(I(:) == i);
    hist_cumul(i) = sum_hist; 
end
hist_max = hist_cumul(255);
map_hist = zeros(1, 256);
j = 0;
seuil_inf = 0;
seuil_sup = hist_max/255;
for i = 1:255
    val = hist_cumul(i);
    while val > seuil_sup
        j = j + 1;
        seuil_inf = seuil_sup;
        seuil_sup = j * hist_max/255;
    end
    (val - seuil_inf) / (seuil_sup - seuil_inf);
    map_hist(i+1) = j  + round((val - seuil_inf) / (seuil_sup - seuil_inf));
end
end