n_games = 165*50;
x = rand(1,n_games)>0.5;

n_winning_streaks = zeros(1,100);
streak = 0;
for i = 1:n_games
    if x(i)
        streak = streak+1;
        n_winning_streaks(streak) = n_winning_streaks(streak)+1;
    else
        streak = 0;
    end
end

plot(1:12, n_winning_streaks(1:12))