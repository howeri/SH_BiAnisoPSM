function plotImage(phi,theta,y,plotName)
imagesc(phi, theta, y)
axis equal
axis tight
xlabel('\theta (degree)')
xticks([0 90 170])
ylabel('\phi (degree)')
yticks([0 90 180 270 350])
set(gca,'FontSize',20)
colorbar
title(plotName)
end

