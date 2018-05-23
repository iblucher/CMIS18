resolutions = [0.1, 0.3, 0.5];
ycoordv1 = [6.841, 2.2937, 1.7328, 1.3259, 1.2859];
ycoordv2 = [1.3301, 1.277, 1.2207];
plot(resolutions, ycoordv2);
xlabel('Mesh resolution');
ylabel('Absolute vaulue of y coordinate');
ylim([0 4]);