FROM torreya/platform_1:latest
WORKDIR /usr/src/

COPY . .

RUN mkdir RLT_files && mkdir RLT_files/full && mkdir RLT_files/simplified && rm -r entry
RUN find ./CodePile -type f -print0 | xargs -0 dos2unix
RUN mkdir /usr/src/WorkSpace && g++ CodePile/rlt_pile/generate.cpp -o /usr/src/CodePile/bin/GRLT -std=c++17 && g++ CodePile/effc_pile/launch.cpp -o /usr/src/CodePile/bin/Geffc -std=c++17
RUN chmod -R 777 /usr/src/CodePile && ln -s /usr/src/CodePile/bin/All_in_ENTRY.py CodePile/bin/MASS && ln -s /usr/src/CodePile/sim_pile/titrate.py CodePile/bin/ODE_solver
ENV PATH="${PATH}:/usr/src/CodePile/bin"

CMD ["/bin/sh"]
