function plot_channel(f, Chann, tit)

    Channel = abs(Chann);
    figure
    grid on
    plot(f, Channel)
    xlabel("Frequency (Hz)")
    ylabel("|Channel|")
    title(tit)
end