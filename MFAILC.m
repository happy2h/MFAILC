% 期望轨迹
for k = 1:1:500
	if k < 250
		yd(k+1) = 0.5*(-1).^(round(k/100));
	else 
		yd(k+1) = 0.5*sin((k*pi)/100) + 0.3*cos((k*pi)/50);
    end
end

% 参数设置
epsilon = 0.01; % 拟伪偏导重置阈值
eta = 1; % 拟伪偏导更新增益
rho = 0.2; % 控制律增益
lamda = 1; % 控制律步长
mu = 2; % 拟伪偏导更新步长

% 控制过程
i_n = 60;  %迭代次数
y(1:i_n,1:500) = 0;
for i = 1:1:i_n
	for k = 1:1:500
		if k == 1
			phi(i,k) = 0.4; % 初始拟伪偏导
		elseif k == 2
			phi(i,k) = phi(i,k-1) + (eta*(ub(i,k-1) - 0)/(mu + norm(ub(i,k-1) - 0)^2))*(y(i,k) - 0 - phi(i,k-1)*(ub(i,k-1) - 0));
		else 
			phi(i,k) = phi(i,k-1) + (eta*(ub(i,k-1) - ub(i,k-2))/(mu + norm(ub(i,k-1) - ub(i,k-2))^2))*(y(i,k) - y(i,k-1) - phi(i,k-1)*(ub(i,k-1) - ub(i,k-2)));
		end
		if i == 1
			uf(i,k) = 0; % 初始控制信号
		else 
			uf(i,k) = uf(i-1,k) + 0.4*e(i-1,k+1);
		end
		if k == 1
			ub(i,k) = 0; % 初始控制信号
		else 
			ub(i,k) = ub(i,k-1) + (rho*phi(i,k)/(lamda + norm(phi(i,k))^2))*(yd(k+1) - y(i,k));
		end
		if k>2 && (phi(i,k) <= epsilon || (abs(ub(i,k-1) - ub(i,k-2)) <= epsilon))
			phi(i,k) = phi(i,1);
		end
		u(i,k) = uf(i,k) + ub(i,k);
		% 系统函数
		if k < 250 
			y(i,k+1) = y(i,k)*u(i,k)/(1 + norm(y(i,k))^2) + (u(i,k) + (1+*round(k/500))*sin(y(i,k)))^3;
		else
			y(i,k+1) = y(i,k)*u(i,k)^3/(1 + norm(y(i,k))^2) + u(i,k)^3;
		end
		e(i,k+1) = yd(k+1) - y(i,k+1);
	end
end
%误差
for i =1:1:i_n
e_min(i) = max(abs(e(i,:)));
end

figure(1) 
plot(yd,'r'); hold on;
plot(y(i_n,:),'b'); title('µü´ú10´Î');
figure(2)
plot(e_min);title('error of time k');
