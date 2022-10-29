from re import L, X
from manim import *
import numpy as np
from math import *
from manim.mobject.geometry import ArrowTriangleFilledTip
import time
class Minkowski(VGroup):
    
    def __init__(
        self, 
        grid_on = False,
        obs_on = False,
        obs_color = None,
        obs_name = None,
        axes_config = None,
        xrange = [-4, 4, 1],
        yrange = [-4, 4, 1], 
        x_length = 5,
        y_length = 5,
        axes_numbers_on = True,
        t = "t", 
        x = "x",
    ):
        super().__init__()
        ac = {
            "numbers_to_include": np.arange(xrange[0], xrange[1]+2, 2),
            "number_scale_value": 0.5,
        }
        if axes_numbers_on == False:
            ac = {}
        if axes_config == None:
            self.grid = Axes(
                x_range =xrange, y_range = yrange, x_length = x_length, y_length = y_length, 
                tips=False,
                axis_config=ac,
                )
        else:
            self.grid = Axes(axes_config)
        y_label = self.grid.get_y_axis_label(t, direction=LEFT, buff = 0.5)
        x_label = self.grid.get_x_axis_label(x)
        self.grid_labels = VGroup(x_label, y_label)
        if obs_on:
            self.line = Line(self.grid.get_y_axis().get_start(), self.grid.get_y_axis().get_end(), color = obs_color)
            self.obs_label = always_redraw(lambda: Tex(f"{obs_name}", color = obs_color).scale(0.8).next_to(self.line.get_top(), RIGHT*0.4))
            self.add(self.grid, self.grid_labels, self.line, self.obs_label)
        else:
            self.add(self.grid, self.grid_labels)
        
        
    def get_grid(self):
        return self.grid
def current_milli_time():#current time in milliseconds
    return round(time.time() * 1000)
class part1(Scene):
    def construct(self):
        t = Tex("Introduction")
        self.play(Create(t))
        self.wait()
class part2(Scene):
    def construct(self):
        t = Tex("The $k$ factor")
        self.play(Create(t))
        self.wait()
class part3(Scene):
    def construct(self):
        t = Tex("Deriving $k$")
        self.play(Create(t))
        self.wait()
class part4(Scene):
    def construct(self):
        t = Tex("Lorentz transformations")
        self.play(Create(t))
        self.wait()
class outro(Scene):
    def construct(self):
        t = Tex(
            r'''That is it for the video. Please keep in mind that this is only the introduction of introductions, 
            and there is so much more to special relativity. \\
            Thanks for watching!'''
            ).scale(0.8)
        self.play(Create(t), run_time = 5)
        self.wait()
class beginning(Scene):
    def construct(self):
        t = Tex(
            '''Special relativity is a very fascinating topic, but very confusing as well. \\
            Today, I will show you a simplified approach to special relativity and use simple algebra to explain this topic.'''
            ).scale(0.8)
        self.play(Create(t), run_time = 4)
        self.wait(2)
class lorentzExplanations(Scene):
    def construct(self):
        lorentz1 = MathTex(r"x' = ",r"\frac{1}{\sqrt{1-v^2}}","(x -vt)", color = BLUE).scale(0.8)
        lorentz2 = MathTex(r"t' = ", r"\frac{1}{\sqrt{1-v^2}}","(t -vx)", color = BLUE).scale(0.8).next_to(lorentz1, DOWN)
        
        lorentzeqs = VGroup(lorentz1, lorentz2)
        box = SurroundingRectangle(lorentzeqs, color=BLUE)
        lorentzeqs += box
        lorentzeqs.move_to([0,0,0])
        gamma1 = MathTex(r"\gamma", color = BLUE).scale(0.8).move_to(lorentz1[1])
        gamma2 = MathTex(r"\gamma", color = BLUE).scale(0.8).move_to(lorentz2[1])
        
        self.add(lorentzeqs)
        self.wait()
        t = Tex(
            r'''
            There we have it, the Lorentz transformations!
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            For simplicity, $\frac{1}{\sqrt{1-v^2}}$ is called $\gamma$, or the Lorentz factor.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.play(Transform(lorentz1[1], gamma1), Transform(lorentz2[1], gamma2))
        self.play(lorentz1[1].animate.shift(LEFT*0.5), lorentz2[1].animate.shift(LEFT*0.5))
        self.play(lorentz1[2].animate.next_to(lorentz1[1], RIGHT), lorentz2[2].animate.next_to(lorentz2[1], RIGHT))
        lorentzeqs = VGroup(lorentz1, lorentz2)
        self.play(Transform(box, SurroundingRectangle(lorentzeqs, color=BLUE)))
        lorentzeqs += box
        self.play(lorentzeqs.animate.move_to([0,0,0]))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            So, what does this all mean? 
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            First of all, we can measure an event and calculate the coordinate of the event measured by someone else moving at constant velocity $v$.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(3)
        projection = MathTex(r"(x,t) \rightarrow (x',t')").shift(DOWN*2)
        self.play(Create(projection))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            In other words, our view of space and time and be transformed into someone else's.\\ 
            Let's see it in action!
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Observer B measures the event at a completely different time and place than observer A does. Also, 
            notice how the speed of light stays the same!
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Many interesting phenomenons occur when we apply the Lorentz transformations. For example, a moving length 
            appears to be contracted to a stationary observer.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Something similar happens with time intervals. If A experiences some time interval, B will measure it to be
            shorter. 
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            In other words, B sees that A's moving clock runs slower. 
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Interestingly, A also sees that B's clock runs slower.  
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            The observations of A and B can seem contradicting, but they are both correct.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            The order of events can also change from different perspectives. 
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        
class lorentz(Scene):
    def construct(self):
        text = Tex("Now, let's use $k$ to derive the Lorentz Transformations.")
        self.play(Create(text))
        self.wait()
        self.play(FadeOut(text))

        mks_a = Minkowski(obs_on = True, obs_color=BLUE, obs_name = "A", xrange=[0, 4, 1], yrange=[0,8,1], x_length=2.3142, 
        y_length= 4.6286, axes_numbers_on=False).shift(DOWN)
        b_function = lambda x: 2*x
        graph_b= mks_a.grid.get_graph( b_function, color=YELLOW, x_range = [0, 4])
        
        b_label = Tex("B", color = YELLOW).next_to(graph_b.get_end(), RIGHT*0.4).scale(0.8)
        light_graph = lambda x: x+1
        light_graph2 = lambda x: -x+7
        p1 = intersection(light_graph, b_function)
        p2 = intersection(light_graph2, b_function)
        light = DashedLine(mks_a.grid.c2p(0,1,0), mks_a.grid.c2p(3,4,0), color = ORANGE)
        light2 = DashedLine(mks_a.grid.c2p(0,7,0), mks_a.grid.c2p(3,4,0), color = ORANGE)
        itst = Dot(point=mks_a.grid.c2p(3,4,0))
        self.play(FadeIn(mks_a, graph_b, b_label, light, light2, itst))
        
        mks_b = Minkowski(obs_on = True, obs_color=YELLOW, obs_name = "B", xrange=[-4,2, 1], yrange=[0,8,1], x_length=3.4713, 
        y_length= 4.6286, axes_numbers_on=False).shift(DOWN).to_edge(buff = 0.5)
        a_function = lambda x: -2*x
        graph_a= mks_b.grid.get_graph( a_function, color=BLUE, x_range = [-4, 0])
        
        a_label = Tex("A", color = BLUE).next_to(graph_a.get_start(), RIGHT*0.4).scale(0.8)
        event = Dot(point=mks_b.grid.c2p(1.1547,2.88675,0))
        light_graph3 = lambda x: x-1.1547 + 2.88675
        light_graph4 = lambda x: -(x-1.1547) + 2.88675
        p1 = intersection(light_graph3, a_function)
        p2 = intersection(light_graph4, a_function)
        light3 = DashedLine(mks_b.grid.c2p(1.1547,2.88675,0), mks_b.grid.c2p(*p1), color = ORANGE)
        light4 = DashedLine(mks_b.grid.c2p(1.1547,2.88675,0), mks_b.grid.c2p(*p2), color = ORANGE)
        self.play(FadeIn(mks_b, graph_a, a_label, light3, light4, event))
        self.wait()
        #draw lines, (x,t), "how do we get from x,t to x',t'? ", page of algebra
        t = Tex(
            r'''
            A and B are trying to measure the position and time of the event represented by the dot.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(3)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            A measures the event to be at $(x, t)$, and for B at $(x', t')$
            '''
             ).shift(UP*2.5).scale(0.7)
        loc_a = Tex("$(x,t)$").next_to(itst, RIGHT+UP, buff = 0).scale(0.7)
        loc_b = Tex("$(x',t')$").next_to(event, RIGHT, buff = 0).scale(0.7)
        self.play(Create(t), Create(loc_a), Create(loc_b))
        self.wait(3)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Given $(x,t)$, how do we find $(x',t')$? This is what the Lorentz Transformations answer. 
            They transform the coordinates $(x,t)$ into $(x',t')$ and enables us to go from A's perspective to B's.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t), run_time = 3)
        self.wait(5)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Using the fact that the speed of light = 1, we can label some times. 
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        t_ = Tex("$t$").next_to(mks_a.grid.c2p(0,4,0), LEFT*0.5).scale(0.7)
        t_minus_x = Tex("$t-x$").next_to(mks_a.grid.c2p(0,1,0), LEFT*0.5, buff = 0).scale(0.7)
        t_plus_x = Tex("$t+x$").next_to(mks_a.grid.c2p(0,7,0), LEFT*0.5, buff = 0).scale(0.7)
        t_prime = Tex("$t'$").next_to(mks_b.grid.c2p(0,2.88675,0), LEFT*0.5).scale(0.7)
        t_prime_plus = Tex("$t'+x'$").next_to(mks_b.grid.c2p(0,4.04145,0), RIGHT*0.5, buff = 0).scale(0.7)
        t_prime_minus = Tex("$t'-x'$").next_to(mks_b.grid.c2p(0,1.73205,0), RIGHT*0.5, buff = 0).scale(0.7)
        self.play(Create(t_), Create(t_minus_x), Create(t_plus_x))
        self.play(Create(t_prime), Create(t_prime_minus), Create(t_prime_plus))
        self.wait(5)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Like when deriving $k$, we can use emission and reception time intervals to write equations.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(3)
        k = Tex(r"$\Delta t_{reception}$", "$=k$", r"$\Delta t_{emission}$").move_to([4.5, 0, 0]).scale(0.8)
        k2 = k.copy()
        self.play(Write(k))
        left = Tex("$t'-x'$").next_to(k[1], LEFT, buff = 0).scale(0.8)
        right = Tex("$(t-x)$").next_to(k[1], RIGHT, buff = 0).scale(0.8)
        self.play(FadeOut(k[0], k[2]),
         ReplacementTransform(t_prime_minus.copy(), left), ReplacementTransform(t_minus_x.copy(), right),
         run_time = 3)
        eq1 = VGroup(left, k[1], right)
        self.wait(3)
        self.play(eq1.animate.shift(UP))
        self.play(Write(k2))
        left2 = Tex("$t+x$").next_to(k2[1], LEFT, buff = 0).scale(0.8)
        right2 = Tex("$(t'+x')$").next_to(k2[1], RIGHT, buff = 0).scale(0.8)
        self.play(FadeOut(k2[0], k2[2]),
         ReplacementTransform(t_plus_x.copy(), left2), ReplacementTransform(t_prime_plus.copy(), right2),
         run_time = 3)
        eq2 = VGroup(left2, k2[1], right2)
        self.wait(3)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Now, time for some "boring algebra" (The steps are simplified).
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()
        self.play(FadeOut(t, mks_a, graph_b, b_label, light, light2, itst,
        mks_b, graph_a, a_label, light3, light4, event, loc_a, loc_b, t_, t_minus_x, t_plus_x, t_prime, 
        t_prime_minus, t_prime_plus))
        self.play(eq1.animate.move_to([0, 3.4, 0]), eq2.animate.move_to([0, 2.8, 0]))
        eq1_ = MathTex(r"kt'-kx' = k^2t-k^2x").scale(0.8).move_to([-4,2,0])
        eq2_ = MathTex(r"kt'+kx' = t+x").scale(0.8).next_to(eq1_, DOWN)
        eq3 = MathTex(r"2kt' =k^2t-k^2x+t+x").scale(0.8).next_to(eq2_, DOWN)
        eq4 = MathTex(r"t' = \frac{k^2+1}{2k}t - \frac{k^2-1}{2k}x").scale(0.8).next_to(eq3, DOWN)
        right_eqs = MathTex(r'''
            
            kt'-kx' &= k^2t-k^2x \\
            kt'+kx' &= t+x \\
            2kt' &=k^2t-k^2x+t+x \\''', 
            r'''
            t' &= \frac{k^2+1}{2k}t - \frac{k^2-1}{2k}x 
            
        ''').scale(0.8).shift(RIGHT*3, UP*0.5)
        eq5 = MathTex(r"2kx' = t+x - k^2t +k^2x").scale(0.8)
        eq6 = MathTex(r"x' = \frac{k^2+1}{2k}x - \frac{k^2-1}{2k}t").scale(0.8)
        left_eqs = MathTex(r'''

            kx + kt' &= t+x \\
            -kx + kt' &= k^2t-k^2x \\
            2kx' &= t+x - k^2t +k^2x \\''',
            r'''
            x' &= \frac{k^2+1}{2k}x - \frac{k^2-1}{2k}t 
            
        ''').scale(0.8).shift(LEFT*3, UP*0.5)
        self.play(Write(left_eqs), run_time = 6)
        self.wait(2)
        self.play(Write(right_eqs), run_time = 6)
        self.wait(2)
        self.play(FadeOut(left_eqs[0], right_eqs[0]))
        self.play(left_eqs[1].animate.shift(UP*2), right_eqs[1].animate.shift(UP*2))
        k = MathTex(r"k = \sqrt{\frac{1+v}{1-v}}").scale(0.8).shift(LEFT*4)
        k1 = MathTex(r"\frac{k^2+1}{2k} = \frac{1}{\sqrt{1-v^2}}").scale(0.8).next_to(k, DOWN)
        k2 = MathTex(r"\frac{k^2-1}{2k} = \frac{v}{\sqrt{1-v^2}}").scale(0.8).next_to(k1, DOWN)
        k.shift(RIGHT*0.5)
        self.play(Write(k))
        self.play(Write(k1))
        self.play(Write(k2))
        lorentz1 = MathTex(r"x' = \frac{1}{\sqrt{1-v^2}}(x -vt)", color = BLUE).scale(0.8).shift(RIGHT*2., DOWN)
        lorentz2 = MathTex(r"t' = \frac{1}{\sqrt{1-v^2}}(t -vx)", color = BLUE).scale(0.8).next_to(lorentz1, DOWN)
        self.play(Write(lorentz1))
        self.play(Write(lorentz2))
        lorentzeqs = VGroup(lorentz1, lorentz2)
        box = SurroundingRectangle(lorentzeqs, color=BLUE)
        lorentzeqs += box
        self.play(Create(box))
        self.play(FadeOut(k, k1, k2, left_eqs[1], right_eqs[1], eq1, eq2))
        self.play(lorentzeqs.animate.move_to([0,0,0]))
        
class dopplerProof(Scene):
    def construct(self):
        t = Tex(
            r'''
            Consider this scenario of B moving away from A at speed $v$ 
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()
        mks_a = Minkowski(obs_on = True, obs_color=BLUE, obs_name = "A", xrange=[0, 4, 1], yrange=[0,7,1], x_length=2.5714, 
        y_length= 4.5, axes_numbers_on=False).shift(DOWN)
        b_function = lambda x: 2*x
        graph_b= mks_a.grid.get_graph( b_function, color=YELLOW, x_range = [0, 3.5])
        
        b_label = Tex("B", color = YELLOW).next_to(graph_b.get_end(), RIGHT*0.4).scale(0.8)
        light_graph = lambda x: x+2
        light_graph2 = lambda x: -x+6
        p1 = intersection(light_graph, b_function)
        p2 = intersection(light_graph2, b_function)
        light = DashedLine(mks_a.grid.c2p(0,2,0), mks_a.grid.c2p(*p1), color = ORANGE)
        light2 = DashedLine(mks_a.grid.c2p(0,6,0), mks_a.grid.c2p(*p2), color = ORANGE)
        itst = Dot(point=mks_a.grid.c2p(*p1))
        self.play(FadeIn(mks_a))
        self.play(Create(graph_b), Create(b_label))
        self.play(Create(light), Create(light2), FadeIn(itst))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            A emits a signal and B reflects it back upon receiving it.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()
        self.play(FadeOut(t))
        t = Tex(
            r'''
            If B receives the signal at time $t$, then we can calculate the times A emits and receives the signal back using
            the fact that the speed of light = 1 in the diagram.
            '''
             ).shift(UP*2.5).scale(0.7)
        time = Tex("$t$").next_to(mks_a.grid.c2p(0,4,0), LEFT*0.5).scale(0.8)
        self.play(Create(t, run_time = 2), Create(time))
        self.wait(2)
        uptime = Tex("$t+d$").next_to(mks_a.grid.c2p(0,6,0), LEFT*0.5).scale(0.8)
        downtime = Tex("$t-d$").next_to(mks_a.grid.c2p(0,2,0), LEFT*0.5).scale(0.8)
        deltax = Line(mks_a.grid.c2p(0,4,0), mks_a.grid.c2p(2,4,0))
        deltax2 = deltax.copy()
        deltax_ = deltax.copy()
        d_label = Tex("$d$").next_to(deltax, UP*0.5)
        self.add(deltax_)
        self.play(Create(deltax), Create(d_label))
        self.play(Rotate(deltax, PI/2, about_point = deltax.get_start()))
        self.play(Create(uptime))
        self.play(Rotate(deltax2, -PI/2, about_point = deltax2.get_start()))
        self.play(Create(downtime))
        self.wait()

        self.play(FadeOut(t))
        k = r"$k = \frac{\Delta t_{reception}}{\Delta t_{emission}}$"
        t = Tex(
            "Recall that ", k
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        k = t[1]
        self.play(FadeOut(t[0]))
        self.play(k.animate.move_to([4.5, 0, 0]))
        mks_b = Minkowski(obs_on = True, obs_color=YELLOW, obs_name = "B", xrange=[-4, 0, 1], yrange=[0,7,1], x_length=2.5714, 
        y_length= 4.5, axes_numbers_on=False).shift(DOWN).to_edge(buff = 0.5)
        a_function = lambda x: -2*x
        graph_a= mks_b.grid.get_graph( a_function, color=BLUE, x_range = [-3.5, 0])
        
        a_label = Tex("A", color = BLUE).next_to(graph_a.get_start(), RIGHT*0.4).scale(0.8)
        light_graph3 = lambda x: x+3.4641
        light_graph4 = lambda x: -x+3.4641
        p1 = intersection(light_graph3, a_function)
        p2 = intersection(light_graph4, a_function)
        light3 = DashedLine(mks_b.grid.c2p(0,3.4641,0), mks_b.grid.c2p(*p1), color = ORANGE)
        light4 = DashedLine(mks_b.grid.c2p(0,3.4641,0), mks_b.grid.c2p(*p2), color = ORANGE)
        self.play(FadeIn(mks_b, graph_a, a_label, light3, light4))

       
        t = Tex(
            r'''
            If we pretend that A emitted a signal at t = 0 that was instantaneously received by B,
            $t-d$ is $\Delta t_{emission}$, and $k(t-d)$ is $\Delta t_{reception}$
            '''
             ).shift(UP*2.5).scale(0.7)
        origin1 = Dot(point=mks_a.grid.c2p(0,0,0))
        origin2 = Dot(point=mks_b.grid.c2p(0,0,0))
        l1 = Line(mks_a.grid.c2p(0,0,0), mks_a.grid.c2p(0,2,0), color = GREEN)
        l2 = Line(mks_b.grid.c2p(0,0,0), mks_b.grid.c2p(0,3.4641,0), color = GREEN)
        label1 = BraceLabel(l1, r"\Delta t_{emission}", brace_direction=RIGHT, label_scale=0.8)
        label2 = BraceLabel(l2, r"\Delta t_{reception}", brace_direction=LEFT, label_scale=0.8)
        t_b = Tex("$k(t-d)$").next_to(mks_b.grid.c2p(0,3.4641,0), RIGHT*0.5).scale(0.8)
        self.play(Create(t, run_time = 2))
        self.wait(4)
        self.play(Create(origin1), Create(origin2))
        self.play(Create(label1))
        self.play(Create(label2))
        self.play(Create(t_b))
        self.wait(3)
        self.play(FadeOut(t, label1, label2))
        l1 = Line(mks_b.grid.c2p(0,0,0), mks_b.grid.c2p(0,3.4641,0), color = GREEN)
        l2 = Line(mks_a.grid.c2p(0,0,0), mks_a.grid.c2p(0,6,0), color = GREEN)
        label1 = BraceLabel(l1, r"\Delta t_{emission}", brace_direction=LEFT, label_scale=0.8)
        label2 = BraceLabel(l2, r"\Delta t_{reception}", brace_direction=RIGHT, label_scale=0.8)
        t = Tex(
            r'''
            If we pretend that B emitted a signal at t = 0 that was instantaneusly received by A,
            now $k(t-d)$ is $\Delta t_{emission}$, and $k \cdot k(t-d)$ is $\Delta t_{reception}$
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t), run_time = 2)
        self.wait(3)
        self.play(Create(label1))
        self.play(Create(label2))
        t_a = Tex("$k^2(t-d)$").next_to(mks_a.grid.c2p(0,6,0), RIGHT*0.5).scale(0.8)
        self.play(Create(t_a))
        self.wait(3)
        self.play(FadeOut(t, label1, label2))
        t = Tex(
            r'''
            We have an equation!
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        eq1 = Tex(r"$t+d = k^2(t-d)$").scale(0.8).move_to(k)
        self.play(Transform(k, eq1))
        self.wait(2)
        self.play(k.animate.shift(UP*1))
        eq2 = Tex(r"$k = \sqrt{\frac{t+d}{t-d}}$").scale(0.8).move_to(eq1)
        self.play(Create(eq2))
        self.play(FadeOut(t))
        self.play(FadeOut(mks_b, graph_a, a_label, light3, light4, t_b, origin2))
        t = Tex(
            r'''
            $d$ is the position of B at time t. 
            Since B is moving at velocity v, $d = vt$
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t), Wiggle(d_label))
        eq3 = Tex(r"$k = \sqrt{\frac{t+vt}{t-vt}}$").scale(0.8).move_to(eq1)
        eq4 = Tex(r"$k = \sqrt{\frac{t(1+v)}{t(1-v)}}$").scale(0.8).move_to(eq1)
        eq_final = Tex(r"$k = \sqrt{\frac{1+v}{1-v}}$").scale(0.8).move_to(eq1)
        self.wait(3)
        d_final = Tex("$vt$").next_to(deltax_, UP*0.5)
        self.play(FadeOut(k))
        self.play(ReplacementTransform(d_label, d_final))
        self.play(ReplacementTransform(eq2, eq3))
        self.wait(2)
        self.play(ReplacementTransform(eq3, eq4))
        self.wait(2)
        self.play(ReplacementTransform(eq4, eq_final))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Great! We have derived the value of k in terms of v!
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()
class dopplerexample(Scene):
    def construct(self):
        t = Tex(
            r'''
            $k$ is known as the relativistic Doppler factor. Because k works for time intervals, it also applies to 
            the period($T$), or the time 1 oscillation takes for a wave. 
            '''
             ).shift(UP*2.5).scale(0.7)
        k = MathTex(r"k = \frac{T_{received}}{T_{emitted}")
        self.play(Create(t))
        self.wait(5)
        self.play(Create(k))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            Because a wave's wavelength($\lambda$) is porportional to its period, $k$ is also the ratio between emitted and received 
            wavelengths.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(4)
        self.play(Transform(k, MathTex(r"k = \frac{cT_{received}}{cT_{emitted}")))
        self.wait(2)
        self.play(Transform(k, MathTex(r"k = \frac{\lambda_{received}}{\lambda_{emitted}")))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            The color of light depends on its wavelength, so moving objects appear to "shift" their color 
            in the color spectrum.
            '''
             ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()



        


        
        
        
        

        
class intro(Scene):
    #fix: diagram used to represent physical events and movements
    def construct(self):
        t = Tex("This is a Minkowski spacetime diagram.").shift(UP*2.5).scale(0.8)
        self.play(Write(t))
        m = Minkowski().shift(DOWN*1)
        self.play(Create(m), run_time = 3)
        self.play(FadeOut(t), FadeOut(m))
        t = Tex("It is kind of like a typical position vs time graph, but with the axes switched").shift(UP*2.5).scale(0.7)
        
        #make normal grid 
        grid = Axes(
                x_range =[-4, 4, 1], y_range = [-4, 4, 1], x_length = 5, y_length = 5, 
                tips=False,
                axis_config={
                    "numbers_to_include": np.arange(-4, 6, 2),
                    "number_scale_value": 0.5,
                },
                )
        y_label = grid.get_y_axis_label("x", direction=LEFT+DOWN*0.5, buff = 0.5)
        x_label = grid.get_x_axis_label("t")
        grid_labels = VGroup(x_label, y_label)
        sin = grid.get_graph(lambda x: np.sin(x), color=BLUE)
        normal = VGroup(grid, grid_labels, sin).shift(DOWN*1)
        self.wait(1)
        self.play(FadeIn(normal))
        self.play(Create(t, run_time = 2))
        self.wait(1)
        self.play(FadeOut(sin))
        normal.remove(sin)
        self.play(Transform(normal, m))
        self.play(FadeOut(t))
        t = Tex("And, rather than graphing functions, Minkowski diagrams represent physical events and movements").shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t =Tex("A point on this diagram is called an event, since it has a position and time.").shift(UP*2.5).scale(0.7)
        p = Dot(point=grid.c2p(1, 1, 0), color=YELLOW)
        p_label = MathTex("(x, t)").next_to(p, UP*0.5).scale(0.8)
        self.play(Create(t, run_time = 2), Create(p))
        self.play(Create(p_label))
        self.wait()
        self.play(FadeOut(t))
        self.play(FadeOut(p, p_label))
        t = Tex(
            "An object staying still would be represented as a vertical line \n since its position doesn’t change with time."
            ).shift(UP*2.5).scale(0.7)
        line = Line(m.grid.c2p(2, -4, 0), m.grid.c2p(2,4,0), color = YELLOW)
        self.play(Create(t))
        self.play(Create(line))
        graph = grid.get_graph(lambda x: -2*x, color=YELLOW, x_range = [-2, 2])
        self.wait()
        self.play(FadeOut(t))
        t = Tex( "A moving object would have a slope.").shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.play(ReplacementTransform(line, graph))
        self.play(FadeOut(t))
        t = Tex("The observer is always staying still at position x = 0, represented by the vertical axis").shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        
        line = Line(m.grid.c2p(0, -4, 0), grid.c2p(0,4,0), color = BLUE)
        a_label = always_redraw(lambda: Tex("A", color = BLUE).next_to(line.get_top(), RIGHT*0.4))
        b_label = always_redraw(lambda: Tex("B", color = YELLOW).next_to(graph.get_start(), LEFT*0.4))
        self.play(Create(line), Create(a_label), Create(b_label))
        self.wait()

        self.play(FadeOut(t))
        t = Tex("In this scenario, according to observer A, object B is moving with a constant velocity to the left.").shift(UP*2.5).scale(0.7)
        a_dot = Dot(color = BLUE).next_to(line.get_top(), UP*0.4)
        b_dot = Dot(color = YELLOW).move_to([1.25, a_dot.get_y(),0])
        #a_label_ = Tex("A", color = BLUE).next_to(a_dot, DOWN*0.4)
        #b_label_ = always_redraw(lambda: Tex("B", color = YELLOW).next_to(b_dot, DOWN*0.4))
        
        self.play(Create(t))
        self.play(Create(a_dot), Create(b_dot))
        self.wait()
        self.remove(graph)
        graph.reverse_direction()
        print(graph.get_direction())
        self.play(b_dot.animate.shift(LEFT*2.5), Create(graph), rate_func = linear, run_time = 2.5)
        self.wait()
        self.play(FadeOut(a_dot, b_dot, t))
        t = Tex(
            '''One important feature of Minkowski diagrams is that the x axis is scaled in light-seconds, 
            or the distance light travels in a second'''
            ).shift(UP*2.5).scale(0.7)
        length = Line(grid.c2p(0,0,0), grid.c2p(-1,0,0), color = GREEN)
        ls = BraceLabel(length, "299,792,458[m]", label_scale = 0.5)
        self.play(Create(t), run_time = 2)
        self.wait()
        self.play(Create(length))
        self.play(Create(ls))

        self.play(FadeOut(t))
        self.play(FadeOut(length, ls))
        t = Tex(
            ''' As a result, a photon is represented by a line tilted 45 degrees
             and any velocity from the diagram is in units of lightseconds/second (light has a speed of 1)'''
            ).shift(UP*2.5).scale(0.7)
        
        light = DashedLine(grid.c2p(-4,-4,0), grid.c2p(4,4,0), color = ORANGE)
        angle = Angle(grid.get_x_axis(), light)
        angle_label = MathTex("45deg").scale(0.5).move_to([angle.get_x()+ 0.3, angle.get_y()+0.3, 0]) #TODO
        self.play(Create(t), run_time = 3)
        self.wait(1.5)
        self.play(Create(light), Create(angle), Create(angle_label))
        self.wait()
        
        
class constantSpeed(Scene):
    def construct(self):
        t = Tex(
            "An important postulate in special relativity is the constant speed of light for all observers moving with constant velocities. "
            ).scale(0.7)
        self.play(Create(t), run_time = 2)
        self.wait()
        self.play(FadeOut(t))
        
        t = Tex(
            "Let's look at an example. Here, observers A and B, moving relative to each other, both measure the speed of object C."
            ).scale(0.7)
        self.play(Create(t), run_time = 2)
        mks_a = Minkowski(obs_on = True, obs_color=BLUE, obs_name = "A").shift(DOWN*1.5, LEFT*2.5).scale(0.8)
        mks_b = Minkowski(obs_on = True, obs_color=YELLOW, obs_name = "B").shift(DOWN*1.5, RIGHT*2.5).scale(0.8)
        graph_c_a= mks_a.grid.get_graph(lambda x: 2*x, color=GREEN, x_range = [-2, 2])
        graph_c_b= mks_b.grid.get_graph(lambda x: 3.5*x, color=GREEN, x_range = [-8/7, 8/7])
        
        a_dot = Dot(color = BLUE)
        b_dot = Dot(color = YELLOW).move_to([a_dot.get_x() - 1.5, a_dot.get_y(),0])
        c_dot = Dot(color = GREEN).move_to([a_dot.get_x() - 3, a_dot.get_y(),0])
        a_label = always_redraw(lambda: Tex("A", color = BLUE).scale(0.5).next_to(a_dot, DOWN*0.4))
        b_label = always_redraw(lambda: Tex("B", color = YELLOW).scale(0.5).next_to(b_dot, DOWN*0.4))
        c_label = always_redraw(lambda: Tex("C", color = GREEN).scale(0.5).next_to(c_dot, DOWN*0.4))
        sn = VGroup(a_dot, b_dot, c_dot, a_label, b_label, c_label)
        #sn.save_state()
        self.wait()
        self.play(FadeOut(t))
        self.play(Create(sn))
        self.wait()
        self.play(b_dot.animate.shift(RIGHT*3), c_dot.animate.shift(RIGHT*6), run_time = 3, rate_func = linear)
        self.wait()
        self.play(FadeOut(sn))
        #A's pov
        a_dot = Dot(color = BLUE).move_to([mks_a.grid.c2p(0,0,0)[0],  1,0])
        b_dot = Dot(color = YELLOW).move_to([mks_a.grid.c2p(-1,0,0)[0], a_dot.get_y(),0])
        c_dot = Dot(color = GREEN).move_to([mks_a.grid.c2p(-2,0,0)[0], a_dot.get_y(),0])
        a_label = always_redraw(lambda: Tex("A", color = BLUE).scale(0.5).next_to(a_dot, DOWN*0.4))
        b_label = always_redraw(lambda: Tex("B", color = YELLOW).scale(0.5).next_to(b_dot, DOWN*0.4))
        c_label = always_redraw(lambda: Tex("C", color = GREEN).scale(0.5).next_to(c_dot, DOWN*0.4))
        b_dot.save_state(), b_label.save_state()
        sn = VGroup(a_dot, b_dot, c_dot, a_label, b_label, c_label)
        #B's pov
        b_dot_ = Dot(color = YELLOW).move_to([mks_b.grid.c2p(0,0,0)[0], 1,0])
        a_dot_ = Dot(color = BLUE).move_to([mks_b.grid.c2p(1,0,0)[0], b_dot_.get_y(),0])
        c_dot_ = Dot(color = GREEN).move_to([mks_b.grid.c2p(-8/7,0,0)[0], b_dot_.get_y(),0])
        a_label_ = always_redraw(lambda: Tex("A", color = BLUE).scale(0.5).next_to(a_dot_, DOWN*0.4))
        b_label_ = always_redraw(lambda: Tex("B", color = YELLOW).scale(0.5).next_to(b_dot_, DOWN*0.4))
        c_label_ = always_redraw(lambda: Tex("C", color = GREEN).scale(0.5).next_to(c_dot_, DOWN*0.4))
        a_dot_.save_state(), a_label_.save_state()
        sn2 = VGroup(a_dot_, b_dot_, c_dot_, a_label_, b_label_, c_label_)
        self.play(Create(sn), Create(sn2))
        #m=0.2857
        label_c_a = Tex("C", color = GREEN).next_to(graph_c_a.get_end(), RIGHT*0.4).scale(0.8)
        label_c_b = Tex("C", color = GREEN).next_to(graph_c_b.get_end(), RIGHT*0.4).scale(0.8)
        

        self.play(Create(mks_a), Create(mks_b), run_time = 3, lag_ratio = 0.1 )
        #self.play( run_time = 2)
        #sn.restore()
        self.wait()
        self.play(
            Create(graph_c_a), b_dot.animate.move_to([mks_a.grid.c2p(1,0,0)[0], a_dot.get_y(),0]), 
            c_dot.animate.move_to([mks_a.grid.c2p(2,0,0)[0], a_dot.get_y(),0]), run_time = 3, rate_func = linear
        )
        self.play( 
            Create(graph_c_b), a_dot_.animate.move_to([mks_b.grid.c2p(-1,0,0)[0], b_dot_.get_y(),0]), 
            c_dot_.animate.move_to([mks_b.grid.c2p(8/7,0,0)[0], a_dot.get_y(),0]), run_time = 3, rate_func = linear
        )
        
        self.play(Create(label_c_a), Create(label_c_b))
        
        t = Tex(
            "In this scenario, the speed of object C appears to be different to observers A and B"
            ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            "B measures C's speed to be slower (steeper in the Minkowski diagram) than A's measurement because B is moving relative to A."
            ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait(2)
        self.play(FadeOut(t))
        t = Tex(
            '''However, if object C is a photon, observers A and B, or anyone moving with 
            a constant velocity, will measure its speed to be the constant c, or 1 in the diagram.
            '''
            ).shift(UP*2.5).scale(0.6)
        self.play(Create(t), run_time = 2)
        self.play(FadeOut(graph_c_a, graph_c_b, label_c_a, label_c_b, c_dot, c_label, c_dot_, c_label_))

        light = DashedLine(mks_a.grid.c2p(-4,-4,0), mks_a.grid.c2p(4,4,0), color = ORANGE)
        light2 = DashedLine(mks_b.grid.c2p(-4,-4,0), mks_b.grid.c2p(4,4,0), color = ORANGE)
        c_dot = Dot(color = ORANGE).move_to([mks_a.grid.c2p(-4,0,0)[0], a_dot.get_y(),0])
        c_label = always_redraw(lambda: Tex("C", color = ORANGE).scale(0.5).next_to(c_dot, DOWN*0.4))
        c_dot_ = Dot(color = ORANGE).move_to([mks_b.grid.c2p(-4,0,0)[0], b_dot_.get_y(),0])
        c_label_ = always_redraw(lambda: Tex("C", color = ORANGE).scale(0.5).next_to(c_dot_, DOWN*0.4))
        b_dot.restore(), b_label.restore(), a_dot_.restore(), a_label_.restore() 
        self.play(Create(c_dot), Create(c_dot_))
        self.play(Create(c_label), Create(c_label_))
        self.wait(2)
        self.play(
            Create(light),  b_dot.animate.move_to([mks_a.grid.c2p(1,0,0)[0], a_dot.get_y(),0]), 
            c_dot.animate.move_to([mks_a.grid.c2p(4,0,0)[0], a_dot.get_y(),0]), run_time = 3, rate_func = linear
        )
        self.wait()
        self.play( 
            Create(light2), a_dot_.animate.move_to([mks_b.grid.c2p(-1,0,0)[0], b_dot_.get_y(),0]), 
            c_dot_.animate.move_to([mks_b.grid.c2p(4,0,0)[0], a_dot.get_y(),0]), run_time = 3, rate_func = linear
        )
        self.play(FadeOut(t))
        t = Tex(
            '''
            As we can see, light behaves very differently than ordinary objects.
            '''
            ).shift(UP*2.5).scale(0.8)
        self.play(Create(t))
        self.wait()
        
        self.wait(3)
        '''
        tracker = ValueTracker(0)
        self.play(tracker.animate.set_value(5), rate_func = linear)
        '''
def intersection(f1, f2):
    m1 = f1(1) - f1(0)
    b1 = f1(0)
    m2 = f2(1) - f2(0)
    b2 = f2(0)
    xi = (b1-b2) / (m2-m1)
    yi = m1 * xi + b1
    return [xi, yi, 0]
def gamma(v):
    return 1/(sqrt(1-v*v))

class doppler(Scene):
    def construct(self):
        t = Tex(
            r'''
            Suppose we have observer A and object B.
            A emits 2 light signals separated by time interval $\Delta t$
            '''
            , ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()
        mks_a = Minkowski(obs_on = True, obs_color=BLUE, obs_name = "A").shift(DOWN*1.2)
        b_function = lambda x: 2*x
        graph_b= mks_a.grid.get_graph( b_function, color=YELLOW, x_range = [-2, 2])
        
        b_label = Tex("B", color = YELLOW).next_to(graph_b.get_end(), RIGHT*0.4).scale(0.8)
        light_graph = lambda x: x+1
        light_graph2 = lambda x: x+2
        p1 = intersection(light_graph, b_function)
        p2 = intersection(light_graph2, b_function)
        light = DashedLine(mks_a.grid.c2p(0,1,0), mks_a.grid.c2p(*p1), color = ORANGE)
        light2 = DashedLine(mks_a.grid.c2p(0,2,0), mks_a.grid.c2p(*p2), color = ORANGE)
        self.play(FadeIn(mks_a))
        self.play(Create(graph_b), Create(b_label))
        self.play(Create(light), Create(light2))
        delta1 = Line(mks_a.grid.c2p(0,1,0), mks_a.grid.c2p(0,2,0))
        delta2 = Line(mks_a.grid.c2p(2,2,0), mks_a.grid.c2p(2,4,0), stroke_width=0)
        brace1 = Brace(delta1, buff = 0.1, direction = LEFT, sharpness = 2, stroke_width=0.1)
        brace2 = Brace(delta2, buff = 0.1, direction = RIGHT, sharpness = 2, stroke_width=0.1)
        label1 = MathTex(r"\Delta t").next_to(brace1, LEFT*0.3).scale(0.7)
        label2 = MathTex(r"\Delta t_b").next_to(brace2, RIGHT*0.3).scale(0.7)
        group1 = VGroup(mks_a, graph_b, b_label, light, light2, delta1, delta2, brace1, brace2, label1, label2)
        self.play(Create(delta1), Create(brace1), Create(label1))
        self.play( Create(brace2), Create(label2))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            From A's perspective, B receives the signals separated by time $\Delta t_b$.
            '''
            , ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.wait()
        self.play(FadeOut(t))
        t = Tex(
            r'''
            However, from B's perspective, this time interval is different.
            '''
            , ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.play(group1.animate.to_edge(LEFT, buff=0.5))
        mks_b = Minkowski(obs_on = True, obs_color=YELLOW, obs_name = "B").next_to(mks_a, RIGHT)
        a_function = lambda x: -2*x
        graph_a= mks_b.grid.get_graph( a_function, color=BLUE, x_range = [-2, 2])
        
        a_label = Tex("A", color = BLUE).next_to(graph_a.get_start(), RIGHT*0.4).scale(0.8)
        light_graph = lambda x: x+1.73
        light_graph2 = lambda x: x+3.46
        p1 = intersection(light_graph, a_function)
        p2 = intersection(light_graph2, a_function)
        light_ = DashedLine(mks_b.grid.c2p(0,1.73,0), mks_b.grid.c2p(*p1), color = ORANGE)
        light2_ = DashedLine(mks_b.grid.c2p(0,3.46,0), mks_b.grid.c2p(*p2), color = ORANGE)
        self.play(FadeIn(mks_b, graph_a, a_label, light_, light2_))
        delta3 = Line(mks_b.grid.c2p(0,1.73,0), mks_b.grid.c2p(0,3.46,0))
        brace3 = Brace(delta3, buff = 0.1, direction = RIGHT, sharpness = 2, stroke_width=0.1)
        label3 = MathTex(r"\Delta t'").next_to(brace3, RIGHT*0.3).scale(0.7)
        self.play(Create(delta3), Create(brace3), Create(label3))
        self.play(FadeOut(t))
        t = Tex(
            r'''
            So, what is the ratio between $\Delta t_{emission}$ and $\Delta t_{reception}$?
            '''
            , ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        div = MathTex(r"\rule{0.7cm}{0.4pt}").shift(RIGHT*5.3, UP*1.5)
        self.play(FadeIn(div))
        unit_length = mks_a.grid.coords_to_point(1,0,0)[0] - mks_a.grid.coords_to_point(0,0,0)[0]
        l1 = Line().set_length(1*unit_length).next_to(div, DOWN)
        l2 = Line().set_length(1.73*unit_length).next_to(div, UP)
        self.play(ReplacementTransform(delta1, l1), ReplacementTransform(delta3, l2), run_time = 2)
        top = MathTex(r"\Delta t_{e}").move_to(l2)
        bottom = MathTex(r"\Delta t_{r}").move_to(l1)
        self.wait(2)
        self.play(Transform(l2, top), Transform(l1, bottom), run_time = 1.5)
        k = MathTex("k = ").next_to(div, LEFT*0.3)
        self.play(FadeOut(t))
        t = Tex(
            r'''
            We will call it $k$
            '''
            , ).shift(UP*2.5).scale(0.7)
        self.play(Create(t))
        self.play(FadeIn(k))
        self.wait()





        

'''
function convert()
{
    var input = document.getElementById("in").value;
    input = trim1(input)

    var w = parseFloat(input);

    if (w >= 380 && w < 440)
    {
        red   = -(w - 440) / (440 - 380);
        green = 0.0;
        blue  = 1.0;
    }
    else if (w >= 440 && w < 490)
    {
        red   = 0.0;
        green = (w - 440) / (490 - 440);
        blue  = 1.0;
    }
    else if (w >= 490 && w < 510)
    {
        red   = 0.0;
        green = 1.0;
        blue  = -(w - 510) / (510 - 490);
    }
    else if (w >= 510 && w < 580)
    {
        red   = (w - 510) / (580 - 510);
        green = 1.0;
        blue  = 0.0;
    }
    else if (w >= 580 && w < 645)
    {
        red   = 1.0;
        green = -(w - 645) / (645 - 580);
        blue  = 0.0;
    }
    else if (w >= 645 && w < 781)
    {
        red   = 1.0;
        green = 0.0;
        blue  = 0.0;
    }
    else
    {
        red   = 0.0;
        green = 0.0;
        blue  = 0.0;
    }
    

    // Let the intensity fall off near the vision limits

    if (w >= 380 && w < 420)
        factor = 0.3 + 0.7*(w - 380) / (420 - 380);
    else if (w >= 420 && w < 701)
        factor = 1.0;
    else if (w >= 701 && w < 781)
        factor = 0.3 + 0.7*(780 - w) / (780 - 700);
    else
        factor = 0.0;

    var gamma = 0.80;
    var R = (red   > 0 ? 255*Math.pow(red   * factor, gamma) : 0);
    var G = (green > 0 ? 255*Math.pow(green * factor, gamma) : 0);
    var B = (blue  > 0 ? 255*Math.pow(blue  * factor, gamma) : 0); 

    var hex = "#" + decimalToHex(R) + decimalToHex(G) + decimalToHex(B);
    var output = "RGB value: " + hex;
    document.getElementById("result").firstChild.nodeValue = output;
    document.getElementById("result").style.color = hex;
}

// -->

        
'''

class example1(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=True,
        )
       
    
    def construct(self):
        v = 0.5
        g = 1/(sqrt(1-v*v))
        matrix = np.dot(g, [[1, -v], [-v, 1]])


        def f2(x):
            return 1/v*x
        #graph= self.plane.get_graph(f,x_min=0,x_max=3, color=YELLOW)
        graph2 = self.plane.get_graph(f2, color = YELLOW)
        
        t = Tex("t(s)").shift(UP*3.5, LEFT*0.5).scale(0.7)
      
        x = Tex("x(lightseconds)").shift(RIGHT*5, UP*0.2).scale(0.5)
        self.add(t,x)
        light = DashedLine(ORIGIN, [7,7,0], color = ORANGE)
        light2 = DashedLine(ORIGIN, [-7,7,0], color = ORANGE)
        
        
        '''
        v2 = Line(ORIGIN, [5,2,0], stroke_width= 0)
        
        
        dot2 = always_redraw(lambda: Dot(point = v2.get_end()))'''
        v1 = Line( ORIGIN, [2,1,0], stroke_width= 0)
        dot1 = always_redraw(lambda: Dot(point = v1.get_end()))
        dot_label = always_redraw(lambda: Tex(f"({round(dot1.get_x(), 2)}, {round(dot1.get_y(), 2)})").next_to(v1.get_end(), RIGHT+UP, buff = 0))
        v3 = Line(ORIGIN, UP*2.5, stroke_width= 0)
        v4 = Line(ORIGIN, [1.5, f2(1.5), 0], stroke_width= 0)
        a_line = Line([0,-15,0], [0, 15, 0], color = BLUE)
        label_a = always_redraw(lambda: Tex("A", color = BLUE).next_to(v3.get_end(), LEFT))
        label_b = always_redraw(lambda: Tex("B", color = YELLOW).next_to(v4.get_end(), RIGHT))

        self.add_transformable_mobject(light, light2, graph2, a_line, v1, v3, v4)
        self.add(label_a, label_b, dot1, dot_label)

        
        self.wait(2)
        self.apply_matrix(matrix)


        self.wait()













v = 0.25
g = 1/(sqrt(1-v*v))
class LinearTransformationSceneExample(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=True,
        )
       
    
    def construct(self):


        matrix = np.dot(g, [[1, -v], [-v, 1]])
        self.apply(matrix)

        self.wait()

    def apply(self, matrix):
        self.plane.save_state()

        #self.play(ApplyMethod(basis.shift, DOWN*3))

        def f(x):
            return -1/v*x
        def f2(x):
            return 1/v*x
        #graph= self.plane.get_graph(f,x_min=0,x_max=3, color=YELLOW)
        graph2 = self.plane.get_graph(f2, color = YELLOW)
        graph3 = self.plane.get_graph(lambda x: 2*x , color = GREEN)
        t = Tex("t(s)").shift(UP*3.5, LEFT*0.5).scale(0.7)
      
        x = Tex("x(lightseconds)").shift(RIGHT*5, UP*0.2).scale(0.5)
        self.add(t,x)
        light = DashedLine(ORIGIN, [7,7,0], color = ORANGE)
        light2 = DashedLine(ORIGIN, [-7,7,0], color = ORANGE)
        
        
        v1 = Line( ORIGIN, [2,1,0], stroke_width= 0)
        v2 = Line(ORIGIN, [5,2,0], stroke_width= 0)
        v3 = Line(ORIGIN, UP*2.5, stroke_width= 0)
        v4 = Line(ORIGIN, [1.5, f2(1.5), 0], stroke_width= 0)
        dot1 = always_redraw(lambda: Dot(point = v1.get_end(), color = YELLOW))
        dot2 = always_redraw(lambda: Dot(point = v2.get_end()))
        label_a = always_redraw(lambda: Tex("A").next_to(v3.get_end(), LEFT))
        label_b = always_redraw(lambda: Tex("B", color = YELLOW).next_to(v4.get_end(), RIGHT))
        #object
        obj = Rectangle(width=1.0, height=15, stroke_color = BLUE, fill_color = BLUE, fill_opacity = 0.6).shift(RIGHT*0.5)

        left, right= obj.get_x() -obj.width/2, obj.get_x()+obj.width/2
        left_, right_ = left/g, right/g
        brace_1 = BraceBetweenPoints(RIGHT*left, RIGHT*right)
        l1 = BraceLabel(brace_1, "length = 1")
        brace_2 = BraceBetweenPoints(RIGHT*left_, RIGHT*right_)
        l2 = BraceLabel(Line(start = RIGHT*left_, end = RIGHT*right_), f"L' = {round(right_ - left_, 2)}", label_scale = 0.6)        
        #anim
        group = VGroup(light, light2, v1, v2, v3, v4, graph2, graph3, self.i_hat, self.j_hat, obj)
        self.add_transformable_mobject(group)
        self.add(dot1, dot2, label_a, label_b)
        
        #group.save_state()

        
        #self.play(Create(l1), lag_ratio = 0.1)
        #self.play(FadeOut(l1), lag_ratio = 0.1)        
        
        
        self.apply_matrix(matrix)
        self.play(Create(l2), lag_ratio = 0.1)
        self.wait()
        self.play(FadeOut(l2), lag_ratio = 0.1)

        #self.play(Restore(self.plane), Restore(g),run_time = 3)
        
class ContinuousMotion(Scene):
    
    def construct(self):
        def func(pos):
            return dist(pos, ORIGIN)*(pos[1]*UP+pos[0]*RIGHT)
        
        stream_lines = StreamLines(
            func, 
            x_min=-2, x_max=2, delta_x=0.1,
            y_min=-2, y_max=2, delta_y=0.1,
            padding = 5,
            stroke_width=1,
            virtual_time=2.5,          # use shorter lines
            min_color_scheme_value=1.5,  max_color_scheme_value= 20,
            colors = ['#4245f5', '#4293f5', '#42d7f5', '#42f5bc', '#42f548', '#81f542', '#d7f542', '#f5c242', '#f54242'],
            max_anchors_per_line=80
        )
        self.add(stream_lines)
        stream_lines.start_animation(warm_up=False, flow_speed=1.5)
        self.wait(5)
        #self.play(stream_lines.end_animation())
        #self.wait(stream_lines.virtual_time / stream_lines.flow_speed)

class example2(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=True,
        )
       
    
    def construct(self):
        v = 0.5
        g = 1/(sqrt(1-v*v))
        matrix = np.dot(g, [[1, -v], [-v, 1]])


        def f2(x):
            return 1/v*x
        #graph= self.plane.get_graph(f,x_min=0,x_max=3, color=YELLOW)
        graph2 = self.plane.get_graph(f2, color = YELLOW)
        
        t = Tex("t(s)").shift(UP*3.5, LEFT*0.5).scale(0.7)
      
        x = Tex("x(lightseconds)").shift(RIGHT*5, UP*0.2).scale(0.5)
        self.add(t,x)

        v1 = Line( ORIGIN, [2,1,0], stroke_width= 0)
        v4 = Line(ORIGIN, [1.5, f2(1.5), 0], stroke_width= 0)
        a_line = Line([0,-15,0], [0, 15, 0], color = BLUE)
        v3 = Line(ORIGIN, UP*2.5, stroke_width= 0)
        label_a = always_redraw(lambda: Tex("A", color = BLUE).next_to(v3.get_end(), LEFT))
        label_b = always_redraw(lambda: Tex("B", color = YELLOW).next_to(v4.get_end(), RIGHT))

        obj = Rectangle(width=1.0, height=15, stroke_color = BLUE, fill_color = BLUE, fill_opacity = 0.6).shift(RIGHT*1.5)

        left, right= obj.get_x() -obj.width/2, obj.get_x()+obj.width/2
        left_, right_ = left/g, right/g
        velocity = MathTex("v_b = 0.5").to_edge(UL).add_background_rectangle()
        l1 = BraceLabel(Line(start = RIGHT*left, end = RIGHT*right), "length = 1", label_scale = 0.8)
        l2 = BraceLabel(Line(start = RIGHT*left_, end = RIGHT*right_), f"moving length = {round(right_ - left_, 2)}", brace_direction=UP ,
        label_constructor = Tex, label_scale = 0.8).add_background_rectangle()
        self.add_transformable_mobject( graph2, a_line, v1, v3, v4, obj)
        self.add(label_a, label_b, l1, velocity, l1)

        
        self.wait(2)
        self.remove(l1)
        self.apply_matrix(matrix)


        self.wait()
        self.play(Create(l2))
        self.wait(2)
class example3(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=True,
        )
       
    
    def construct(self):
        v = 0.5
        g = 1/(sqrt(1-v*v))
        matrix = np.dot(g, [[1, -v], [-v, 1]])


        def f2(x):
            return 1/v*x
        #graph= self.plane.get_graph(f,x_min=0,x_max=3, color=YELLOW)
        graph2 = self.plane.get_graph(f2, color = YELLOW)
        
        t = Tex("t(s)").shift(UP*3.5, LEFT*0.5).scale(0.7)
      
        x = Tex("x(lightseconds)").shift(RIGHT*5, UP*0.2).scale(0.5)
        self.add(t,x)

        v1 = Line( ORIGIN, [2,1,0], stroke_width= 0)
        v4 = Line(ORIGIN, [1.5, f2(1.5), 0], stroke_width= 0)
        a_line = Line([0,-15,0], [0, 15, 0], color = BLUE)
        v3 = Line(ORIGIN, UP*2.5, stroke_width= 0)
        label_a = always_redraw(lambda: Tex("A", color = BLUE).next_to(v3.get_end(), LEFT))
        label_b = always_redraw(lambda: Tex("B", color = YELLOW).next_to(v4.get_end(), RIGHT))

        obj = Rectangle(width=15, height=2, stroke_color = BLUE, fill_color = BLUE, fill_opacity = 0.6).shift(UP*2)

        down, up= obj.get_y() -obj.height/2, obj.get_y()+obj.height/2
        down_, up_ = down/g, up/g
        velocity = MathTex("v_b = 0.5").to_edge(UL).add_background_rectangle()
        l1 = BraceLabel(Line(start = UP*down, end = UP*up), "\Delta t = 2", label_scale = 0.8, brace_direction=RIGHT )
        l2 = BraceLabel(Line(start = UP*down_, end = UP*up_), rf"moving $\Delta t= {round(up_ - down_, 2)}$", brace_direction=RIGHT ,
        label_constructor = Tex, label_scale = 0.8)
        self.add_transformable_mobject( graph2, a_line, v1, v3, v4, obj)
        self.add(label_a, label_b, l1, velocity, l1)

        
        self.wait(2)
        self.remove(l1)
        self.apply_matrix(matrix)


        self.wait()
        self.play(Create(l2))
        self.wait(2)
class example4(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=True,
        )
       
    
    def construct(self):
        v = -0.5
        g = 1/(sqrt(1-v*v))
        matrix = np.dot(g, [[1, -v], [-v, 1]])


        def f2(x):
            return 1/v*x
        #graph= self.plane.get_graph(f,x_min=0,x_max=3, color=YELLOW)
        graph2 = self.plane.get_graph(f2, color = BLUE)
        
        t = Tex("t(s)").shift(UP*3.5, LEFT*0.5).scale(0.7)
      
        x = Tex("x(lightseconds)").shift(RIGHT*5, UP*0.2).scale(0.5)
        self.add(t,x)

        v1 = Line( ORIGIN, [2,1,0], stroke_width= 0)
        v4 = Line(ORIGIN, [-1.5, f2(-1.5), 0], stroke_width= 0)
        b_line = Line([0,-15,0], [0, 15, 0], color = YELLOW)
        v3 = Line(ORIGIN, UP*2.5, stroke_width= 0)
        label_a = always_redraw(lambda: Tex("A", color = BLUE).next_to(v4.get_end(), LEFT))
        label_b = always_redraw(lambda: Tex("B", color = YELLOW).next_to(v3.get_end(), RIGHT))

        obj = Rectangle(width=15, height=2, stroke_color = YELLOW, fill_color = YELLOW, fill_opacity = 0.4).shift(UP*2)

        down, up= obj.get_y() -obj.height/2, obj.get_y()+obj.height/2
        down_, up_ = down/g, up/g
        velocity = MathTex("v_a = -0.5").to_edge(UL).add_background_rectangle()
        l1 = BraceLabel(Line(start = UP*down, end = UP*up), "\Delta t = 2", label_scale = 0.8, brace_direction=RIGHT )
        l2 = BraceLabel(Line(start = UP*down_, end = UP*up_), rf"moving $\Delta t= {round(up_ - down_, 2)}$", brace_direction=RIGHT ,
        label_constructor = Tex, label_scale = 0.8)
        self.add_transformable_mobject( graph2, b_line, v1, v3, v4, obj)
        self.add(label_a, label_b, l1, velocity, l1)

        
        self.wait(2)
        self.remove(l1)
        self.apply_matrix(matrix)


        self.wait()
        self.play(Create(l2))
        self.wait(2)
class example5(LinearTransformationScene):
    def __init__(self):
        LinearTransformationScene.__init__(
            self,
            show_coordinates=True,
            leave_ghost_vectors=True,
        )
       
    
    def construct(self):
        v = -0.5
        g = 1/(sqrt(1-v*v))
        matrix = np.dot(g, [[1, -v], [-v, 1]])


        def f2(x):
            return 1/v*x
        #graph= self.plane.get_graph(f,x_min=0,x_max=3, color=YELLOW)
        graph2 = self.plane.get_graph(f2, color = BLUE)
        
        t = Tex("t(s)").shift(UP*3.5, LEFT*0.5).scale(0.7)
      
        x = Tex("x(lightseconds)").shift(RIGHT*5, UP*0.2).scale(0.5)
        self.add(t,x)
        
        v1 = Line( ORIGIN, [-2,1,0], stroke_width= 0)
        v2 = Line(ORIGIN, [-5,2,0], stroke_width= 0)
        dot1 = always_redraw(lambda: Dot(point = v1.get_end(), color = GREEN))
        dot2 = always_redraw(lambda: Dot(point = v2.get_end(), color = RED))
        def get_dot1():
            if dot1.get_y() < dot2.get_y():
                return "earlier"
            else: 
                return "later"
        def get_dot2():
            if dot1.get_y() < dot2.get_y():
                return "later"
            else: 
                return "earlier"
        dot_label = always_redraw(lambda: Tex(f"{get_dot1()}", color = GREEN).next_to(v1.get_end(), RIGHT+UP, buff = 0).scale(0.8))
        dot_label2 = always_redraw(lambda: Tex(f"{get_dot2()}", color = RED).next_to(v2.get_end(), RIGHT+UP, buff = 0).scale(0.8))
        v4 = Line(ORIGIN, [-1.5, f2(-1.5), 0], stroke_width= 0)
        b_line = Line([0,-15,0], [0, 15, 0], color = YELLOW)
        v3 = Line(ORIGIN, UP*2.5, stroke_width= 0)
        label_a = always_redraw(lambda: Tex("A", color = BLUE).next_to(v4.get_end(), LEFT))
        label_b = always_redraw(lambda: Tex("B", color = YELLOW).next_to(v3.get_end(), RIGHT))

        velocity = MathTex("v_a = -0.5").to_edge(UL).add_background_rectangle()

        self.add_transformable_mobject( graph2, b_line, v1, v2, v3, v4)
        self.add(label_a, label_b, velocity, dot_label, dot_label2, dot1, dot2)

        
        self.wait(2)
        self.apply_matrix(matrix)


        self.wait()
        self.wait(2)
        