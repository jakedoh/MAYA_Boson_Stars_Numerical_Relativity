#ifndef FUNCTORS_MONAD_HDR
#define FUNCTORS_MONAD_HDR
#include <functional>
#include <type_traits>
#include <vector>
#include <numeric>
#include <algorithm>

//Note: EitherT, MaybeT, WriterT not done yet
	template<typename F, typename A>
constexpr auto operator>>(A a, F f) -> decltype(f(a))
{
	return f(a);
}
	
	template<typename F, typename A>
constexpr auto operator<<(F f, A a) -> decltype(f(a))
{
	return f(a);
}

template<typename F, typename G>
struct ComposeFunctor
{
	F first;
	G second;
	constexpr ComposeFunctor(F f, G g) : first(f), second(g) {}

	template<typename A>
	constexpr auto operator()(A a) const -> decltype(second(first(a)))
	{
		return second(first(a));
	}
};

//doodle
template<typename F, typename G>
constexpr ComposeFunctor<F,G> operator*(F f, G g)
{
	return ComposeFunctor<F,G>(f, g);
}

namespace Functional
{
	struct Nothing {};

	
	template<typename F, typename A>
	class CurryFunctorHelper2
	{
		private:
			F func;
			A arg1;

		public:
			constexpr CurryFunctorHelper2(F f, A a) : func(f), arg1(a) {}
			constexpr CurryFunctorHelper2() : func(F()), arg1(A()) {}

			template<typename... Args>
			constexpr auto operator()(Args... args) const -> decltype(func(arg1, args...)) 
			{
				return func(arg1, args...);
			}

			/*
		template<typename... Args>
				struct result
				{
					typedef typename std::result_of<F(A, Args...)>::type type;
				};
				*/
	};



	template<typename F>
	class CurryFunctorHelper1
	{
		private:
			F func;

		public:
		constexpr CurryFunctorHelper1(const F& f) : func(f) {}
		constexpr CurryFunctorHelper1() : func(F()) {}

		template<typename Arg1>
			constexpr CurryFunctorHelper2<F,Arg1> operator()(Arg1 a) const
			{
				return CurryFunctorHelper2<F,Arg1>(func, a);
			}

		template<typename Arg1>
				struct result
				{
					typedef CurryFunctorHelper2<F,Arg1> type;
				};
	};

	struct CurryFunctor
	{
		template<typename F>
			constexpr CurryFunctorHelper1<F> operator()(const F& func) const
			{
				return CurryFunctorHelper1<F>(func);
			}
	};

	constexpr CurryFunctor curry = CurryFunctor();
	
	template<template<typename> class Md>
		struct UnitFunctor
		{
			template<typename A>
				constexpr Md<A> operator()(A a) const
				{
					return Md<A>(a);
				}
		};

	template<template<typename> class Md, typename F>
		struct FmapFunctorHelper
		{
			F f;
			constexpr FmapFunctorHelper(const F& fc) : f(fc) {}
			constexpr FmapFunctorHelper() : f(F()) {}

			template<typename A>
				constexpr auto operator()(const Md<A>& ma) const -> Md<decltype(f(A()))>
				{
					return Md<A>::fmap(f, ma);
				}
		};

	template<template<typename> class Md>
		struct FmapFunctor
		{
			template<typename F>
				constexpr FmapFunctorHelper<Md, F> operator()(const F& f) const
				{
					return FmapFunctorHelper<Md, F>(f);
				}
		};


	template<template<typename> class Md>
		struct JoinFunctor
		{
			template<typename A>
				constexpr Md<A> operator()(const Md<Md<A> >& ma) const
				{
					return Md<A>::join(ma);
				}
		};

	template<template<typename> class Md, typename F>
	class BindFunctorHelper
	{
		private:
			F func;

		public:
		constexpr BindFunctorHelper(const F& f) : func(f) {}
		constexpr BindFunctorHelper() : func(F()) {}

		template<typename A>
			constexpr auto operator()(const Md<A>& ma) const -> decltype(func(A())) 
			{
				using B = decltype(func(A()));
				return B::join(Md<A>::fmap(func,ma));	
			}
	};
	
	template<template<typename> class Md>
		struct BindFunctor
		{
			template<typename F>
				constexpr BindFunctorHelper<Md, F> operator()(const F& f) const
				{
					return BindFunctorHelper<Md, F>(f);
				}
		};
	
	template<template<typename> class Md, typename F>
		struct ApplyFunctorHelper
		{
			Md<F> mf;
			constexpr ApplyFunctorHelper(const Md<F>& fc) : mf(fc) {}
			constexpr ApplyFunctorHelper() : mf(Md<F>()) {}

			template<typename A>
				auto operator()(Md<A> ma) const -> Md<decltype(F()(A()))>
				{
					using B = decltype(F()(A()));
					return std::function<Md<B>(F)>([ma](F fop) -> Md<B>
					{
						const auto mapper = FmapFunctor<Md>();
						return mapper << fop << ma;
					}) >> (BindFunctor<Md>()) << mf;
				}
		};

	template<template<typename> class Md>
		struct ApplyFunctor
		{
			constexpr ApplyFunctor() {}
			template<typename F>
				constexpr ApplyFunctorHelper<Md, F> operator()(const Md<F>& mf) const
			{
				return ApplyFunctorHelper<Md, F>(mf);
			}
		};
		
	template<template<typename> class Md, typename B>
		struct ThenFunctorHelper
		{
			Md<B> mb;
			constexpr ThenFunctorHelper(const Md<B>& inb) : mb(inb) {}
			constexpr ThenFunctorHelper() : mb(Md<B>()) {}

			template<typename A>
				Md<B> operator()(Md<A> ma) const
				{
					const auto ourb = mb;
					return std::function<Md<B>(A)>([ourb](A a)
					{
						return ourb;
					}) >> (BindFunctor<Md>()) << ma;
				}
		};

	template<template<typename> class Md>
		struct ThenFunctor
		{
			template<typename B>
				constexpr ThenFunctorHelper<Md, B> operator()(const Md<B>& mb)
				{
					return ThenFunctorHelper<Md,B>(mb);
				}
		};

	
	template<template<typename> class Md>
	struct LiftM2Functor
	{
		//F = A x B -> C
		template<typename F, typename A, typename B>
			constexpr auto operator()(F f, Md<A> a, Md<B> b) const -> Md<decltype(f(A(), B()))> 
			{
/*				const auto opmd = f >> curry >> FmapFunctor<Md>();
				return a >> (b >> (f >> curry >> (FmapFunctor<Md>() >> curry) >> (ApplyFunctor<Md>() >> curry)));
				*/
				return f >> curry >> (UnitFunctor<Md>())  >> (ApplyFunctor<Md>()) << a >> (ApplyFunctor<Md>()) << b;
			}
	};

	template<template<typename> class Md>
		struct LiftM3Functor
		{
			template<typename F, typename A, typename B, typename C>
				constexpr auto operator()(F f, Md<A> a, Md<B> b, Md<C> c) const -> Md<decltype(f(A(), B(), C()))>
				{
					return f >> curry >> FmapFunctor<Md>() << a 
						>> (curry >> FmapFunctor<Md>()) >> ApplyFunctor<Md>() << b
					   	>> ApplyFunctor<Md>() << c;
				}
		};
	
	template<template<typename> class Md>
		struct LiftM4Functor
		{
			template<typename F, typename A, typename B, typename C, typename D>
				constexpr auto operator()(F f, Md<A> a, Md<B> b, Md<C> c, Md<D> d) const -> Md<decltype(f(A(), B(), C(), D()))>
				{
					return f >> curry >> FmapFunctor<Md>() << a 
						>> (curry >> FmapFunctor<Md>()) >> ApplyFunctor<Md>() << b
						>> (curry >> FmapFunctor<Md>()) >> ApplyFunctor<Md>() << c
					   	>> ApplyFunctor<Md>() << d;
				}
		};

	template<template<typename> class Md>
		class Monad
		{
			public:

				static constexpr ApplyFunctor<Md> apply = ApplyFunctor<Md>();
				static constexpr UnitFunctor<Md> unit = UnitFunctor<Md>();
				static constexpr BindFunctor<Md> bind = BindFunctor<Md>();
				static constexpr FmapFunctor<Md> fmap = FmapFunctor<Md>();
				static constexpr JoinFunctor<Md> join = JoinFunctor<Md>();
				static constexpr ThenFunctor<Md> then = ThenFunctor<Md>();
				static constexpr LiftM2Functor<Md> liftM2 = LiftM2Functor<Md>();
				static constexpr LiftM3Functor<Md> liftM3 = LiftM3Functor<Md>();
				static constexpr LiftM4Functor<Md> liftM4 = LiftM4Functor<Md>();
		};

	/*
	template<template<typename> class Md>
	CurryFunctorHelper1<ApplyFunctor<Md> > Monad<Md>::apply;
	*/
	
	template<typename A>
		class Identity
		{
			private:
				A value;

			public:

				constexpr Identity(const A& dt) : value(dt) {}
				constexpr Identity() : value(A()) {}

					constexpr Identity(const Identity<const A>& copy) : value(copy.run()) {}

				constexpr A run() const
				{
					return value;
				}

				constexpr static A run(const Identity<A>& id)
				{
					return id.run();
				}

				constexpr static Identity<A> unit(const A& a)
				{
					return Identity<A>(a);
				}

				template<typename F>
					constexpr static auto fmap(const F& f, const Identity<A>& ia) -> Identity<decltype(f(A()))> 
					{
						return Identity<decltype(f(A()))>(f(ia.run()));
					}

				constexpr static Identity<A> join(const Identity<Identity<A> >& iia)
				{
					return Identity<A>(iia.run().run());
				}
		};


	struct RunIdentityFunctor
	{
		template<typename A>
			constexpr A operator()(const Identity<A>& ia) const
			{
				return ia.run();
			}
	};
	
	constexpr auto runIdentity = RunIdentityFunctor();
	
	template<typename A>
		class Impure
		{
			private:
				std::function<A()> value;

			public:
				Impure(const A& dt)
				{
					value= [dt]() -> A
					{
						return dt;
					};
				}

				Impure(const std::function<A()>& fc)
				{
					value = fc;
				}

				Impure(const Impure<const A>& copy) : value(copy.getValue()) {}

				Impure()
				{
					value = []() { return A(); };
				}

				template<typename F>
					static Impure<A> makeFrom(const F& f)
					{
						return Impure<A>(std::function<A()>(f));
					}


				A run() const
				{
					return value();
				}

				std::function<A()> getValue() const
				{
					return value;
				}

				static A run(const Impure<A>& im)
				{
					return im.run();
				}

				static Impure<A> unit(const A& a)
				{
					return Impure<A>(a);
				}

				template<typename F>
					static auto fmap(const F& f, const Impure<A>& ia) -> Impure<decltype(f(A()))>
					{
						using B = decltype(f(A()));
						const std::function<B()> retf = [f, ia]() -> B
						{
							return f(ia.run());
						};
						return Impure<B>(retf);
					}

				static Impure<A> join(const Impure<Impure<A> >& iia)
				{
					const auto retf = [iia]()
					{
						return iia.run().run();
					};
					return Impure<A>::makeFrom(retf);
				}

		};
		
	template<typename L, typename A>
		class Writer
		{
			private:
				A data;
				L output;

			public:
				constexpr Writer(L out, A val): data(val), output(out) {}
				constexpr Writer(A val): data(val), output(L()) {}
				constexpr Writer() : data(A()), output(L()) {}

					constexpr Writer(const Writer<L, const A>& copy) : data(copy.eval()), output(copy.exec()) {}
				constexpr std::pair<L, A> run() const { return std::make_pair(output, data); } 
				constexpr static std::pair<L, A> run(Writer<L, A> wr) { return wr.run(); }

				constexpr A eval() const { return data; }
				constexpr static A eval(Writer<L,A> wr) { return wr.eval(); }

				constexpr L exec() const { return output; }
				constexpr static L exec(Writer<L,A> wr) { return wr.exec(); }

				template<typename F>
					constexpr static auto fmap(const F& f, const Writer<L,A>& wa) -> Writer<L, decltype(f(A()))>
					{
						return Writer<L, decltype(f(A()))>(wa.output, f(wa.data));
					}

				constexpr static Writer<L,A> join(const Writer<L,Writer<L,A> >& wwa)
				{
					return Writer<L,A>(wwa.exec() + wwa.eval().exec(), wwa.eval().eval());
				}
		};

		struct EvalWriterFunctor
		{
			template<typename W, typename A>
				constexpr A operator()(const Writer<W,A>& wa) const
				{
					return wa.eval();
				}
		};

		constexpr auto evalWriter = EvalWriterFunctor();

		struct ExecWriterFunctor
		{
			template<typename W, typename A>
				constexpr W operator()(const Writer<W,A>& wa) const
				{
					return wa.exec();
				}
		}; 

		constexpr auto execWriter = ExecWriterFunctor();

		struct RunWriterFunctor
		{
			template<typename W, typename A>
				constexpr std::pair<W,A> operator()(const Writer<W,A>& wa) const
				{
					return wa.run();
				}
		};
		
		template<typename L, typename A>
		class Either
		{
			private:
				A right;
				L left;
				bool is_right;
				
				constexpr static Either<L,A> unconst(const Either<L,const A>& copy)
				{
					return copy.is_right ? Either<L,A>(copy.right) : Either<L,A>(copy.left, copy.is_right);
				}

			public:
				constexpr Either(const A& val) : right(val), left(L()), is_right(true) {}
				constexpr Either(const L& l, const bool& correct) : left(l), is_right(correct) {}
				constexpr Either() : right(A()), is_right(true), left(L()) {}
				constexpr Either(const Either<L,const A>& copy) : right(unconst(copy).right), left(unconst(copy).left), is_right(unconst(copy).is_right) {}




				constexpr static Either<L, A> makeRight(const A& val)
				{
					return Either(val);
				}

				constexpr static Either<L, A> makeLeft(const L& le)
				{
					return Either(le, false);	
				}

				constexpr bool isRight() const
				{
					return is_right;
				}

				constexpr static bool isRight(const Either<L,A>& ei)
				{
					return ei.isRight();
				}

				constexpr bool isLeft() const
				{
					return !is_right;
				}

				constexpr static bool isLeft(const Either<L,A>& ei)
				{
					return ei.isLeft();
				}

				template<typename F, typename G>
					constexpr auto either(F leftop, G rightop) const -> decltype(leftop(L()))
					{
						return is_right ? rightop(right) : leftop(left);
					}

				/*
				template<typename B>
					static B either(std::function<B(L)> leftop, std::function<B(A)> rightop, Either<L,A> ei)
					{
						return ei.either(leftop, rightop);
					}
					*/

				template<typename F>
					static auto fmap(const F& f, const Either<L,A>& la) -> Either<L,decltype(f(A()))>
					{
						using B = decltype(f(A()));
							const auto ap = [f](A a) { return Either<L,B>(f(a)); };
							const auto id = [](L l) { return Either<L,B>::makeLeft(l); };
							return la.either(id, ap);
					}

				//this could be made constexpr if, instead of using lambdas, we use some hand-crafted functors.
				static Either<L,A> join(Either<L,Either<L,A> > eea) 
				{
					const auto leftop = [](L l)
					{
						return Either<L,A>::makeLeft(l);
					};

					const auto rightop = [](Either<L,A> ea)
					{
						return ea;
					};

					return eea.either(leftop, rightop);
				}

				bool operator==(Either<L,A> e2)
				{
					const auto cur_is_right = is_right;
					const auto cur_right = right;
					const auto cur_left = left;
					const auto aeq = [cur_is_right, cur_right](A a) { return cur_is_right ? cur_right == a : false; };
					const auto leq = [cur_is_right, cur_left](L l) { return !cur_is_right ? cur_left == l : false; };

					return e2.either(leq, aeq);
				}
		};

		struct IsRightFunctor
		{
			template<typename L, typename A>
				constexpr bool operator()(const Either<L,A>& ea) const
				{
					return ea.isRight();
				};
		};

		constexpr auto isRight = IsRightFunctor();
		
		struct IsLeftFunctor
		{
			template<typename L, typename A>
				bool operator()(const Either<L,A>& ea) const
				{
					return ea.isLeft();
				};
		};
	
	struct EitherFunctor
	{
		template<typename L, typename A, typename F, typename G>
			//notice: the order here is such that currying will require left, then right, then the either. using >> notation, this practically means that it'll look like (data >> (right >> (left >> either))), which may be strange, but using <<, it is (either << left << right << data), or more likely, data >> (either << left << right). This way, the overall dataflow is still left to right, and we can think of either as going left to right as well.
			 constexpr auto operator()(F left, G right, Either<L,A> ea) const -> decltype(left(L()))
			{
				return ea.either(left, right);
			}
	};

	constexpr auto either = (EitherFunctor() >> curry) * curry;
	
	template<typename A>
			class List
			{
				private:
					std::vector<A> values;

				public:
					List() : values(std::vector<A>()) {}
					List(A dt) : values(std::vector<A>(1,dt)) {}
					List(std::vector<A> dts) : values(dts) {}

					std::vector<A> run() const
					{
						return values;
					}


					static std::vector<A> run(List<A> li)
					{
						return li.run();
					}

					static List<A> concat(const List<A>& l1, const List<A>& l2)
					{
						auto d3 = l1.run();
						d3.insert(d3.end(), l2.values.begin(), l2.values.end());
						return List<A>(d3);
					}

					template<typename F>
					static A fold(const F& func, const A& initval, const List<A>& la)
					{
						return std::accumulate(la.values.begin(), la.values.end(), initval, func);
					}

					static List<A> join(const List<List<A> >& ll)
					{
					    return List<List<A> >::fold(*concat, List<A>(), ll);
					}

					int length() const
					{
						return values.size();
					}

					static List<A> range(A start, int len)
					{
						std::vector<A> invec(len);
						std::iota(invec.begin(), invec.end(), start);
						return List<A>(invec);
					}

					A head()
					{
						return values[0];
					}

					template<typename F>
						static auto fmap(const F& f, const List<A>& la) -> List<decltype(f(A()))>
						{
							using B = decltype(f(A()));
							std::vector<B> newvals(la.values.size());
							std::transform(la.values.begin(), la.values.end(), newvals.begin(), f);
							return List<B>(newvals);
						}

					static int length(List<A> ld)
					{
						return ld.length();
					}

					/*
					Maybe<A>& operator[] (const std::size_t idx)
					{
						if(idx >= length() || idx < 0)
						{
							return Maybe<A>();
						}
						else
						{
							return Maybe<A>(values[idx]);
						}
					}
					*/

			};

	/*
	struct RunListFunctor
	{
		template<typename A>
			List<A> operator()(List<A> lta)
			{
				return lta.run();
			}
	};

	constexpr auto runList = RunListFunctor();
	*/
	
	struct HeadListFunctor
	{
		template<typename A>
			A operator()(List<A> lta)
			{
				return lta.head();
			}
	};

	constexpr auto headList = HeadListFunctor();

	struct ConcatFunctor
	{
		template<typename A>
			std::function<List< A>(List<A>)> operator()(List< List<A> > la1) const
			{
				return [la1](List<A> la2)
				{
					return List<A>::concat(la1, la2);
				};
			}
	};
	
	constexpr auto concatList = ConcatFunctor();

	struct FoldListFunctor
	{
		template<typename A, typename F>
//			std::function<std::function<A(List<A>)>(A)> operator()(std::function<A(A,A)> op) const
			A operator()(const F& f, const A& a, const List<A>& la) const
			{
				return List<A>::fold(f, a, la);		
			}
	};
	constexpr auto foldList = (FoldListFunctor() >> curry) * curry;
	
//Stuff for Monad Transformers	
//
	template<template<template<typename> class, typename> class Tr, template<typename> class Md>
		struct LiftFunctor
		{
			template<typename A>
				constexpr Tr<Md, A> operator()(const Md<A>& ma) const
				{
						return Tr<Md, A>::lift(ma);
				}
		};

	template<template<template<typename> class, typename> class Tr, template<typename> class Md>
		class MonadT
		{
			public:
				static constexpr LiftFunctor<Tr, Md> lift = LiftFunctor<Tr, Md>();
		};
	
	template<typename L, template<typename> class Md, typename A>
		class EitherT
		{
			private:
				Md< Either<L,A> > value;
			
			public:
				constexpr EitherT(const A& dt) : value(Md< Either<L, A> >(Either<L, A>(dt))) {}

				constexpr Md<Either<L,A> > run() const
				{
					return value;
				}

				constexpr EitherT(const EitherT<L, Md, const A>& copy) : value(copy.run()) {}

				constexpr Md<Either<L,A> > static run(const EitherT<L, Md, A>& ei)
				{
					return ei.run();
				}

				EitherT(const Md<A>& bd)
				{
					const std::function<Md<Either<L, A> >(A)> bindf =  [](A a) { return EitherT<L,Md,A>(a).run(); };
					value = bd >> (bindf >> Monad<Md>::bind);
				}

				EitherT() : value(Md<Either<L,A> >()) {}

				constexpr static EitherT<L, Md, A> lift(const Md<A>& bd)
				{
					return EitherT<L, Md, A>(bd);
				}

				constexpr static EitherT<L, Md, A> unit(const A& a)
				{
					return EitherT<L, Md, A>(a);
				}

				constexpr EitherT(const Either<L, A>& ei) : value(Md<Either<L,A> >(ei)) {}

				constexpr EitherT(const Md<Either<L, A> >& bme) : value(bme) {}

				constexpr static EitherT<L, Md, A> makeLeft(const L& l)
				{
					return EitherT<L, Md, A>(Either<L,A>::makeLeft(l));
				}

				constexpr static EitherT<L, Md, A> makeRight(const A& a)
				{
					return EitherT<L, Md, A>(Either<L,A>::make_right(a));
				}


				template<typename F>
				static auto fmap(const F& f, const EitherT<L,Md, A>& ea) -> EitherT<L, Md, decltype(f(A())) >
				{
					using B = decltype(f(A()));
					struct LeftOp {
						Either<L,B> operator()(const L& l) const 
						{
							return Either<L, B>::makeLeft(l);
						};
					} leftop;

					struct RightOp {
						F func;
						RightOp(const F& fc) : func(fc) {}

						Either<L,B> operator()(const A& a) const
						{
							return Either<L, B>(func(a));
						}
					} rightop(f);

					return EitherT<L, Md, B>((either << leftop << rightop) >> Monad<Md>::fmap << ea.value);
				}

				static EitherT<L, Md, A> join(const EitherT<L, Md, EitherT<L, Md, A> >& eea);
				/* Okay, our hammer is either: E^L_A x (L->B) x (A->B) -> B
				 * Here, we have M_{E^L_A}
				 * Or rather, M_{EM_A} that we need to turn into M_{E_A}
				 *
				 */
					
		};
	
	struct RunEitherTFunctor
	{
		template<typename L, template<typename> class Md, typename A>
		constexpr Md<Either<L,A> > operator()(EitherT<L, Md, A> ei) const
		{
			return ei.run();
		}
	};

	constexpr auto runEitherT = RunEitherTFunctor();
	
	template<typename L, template<typename> class Md, typename A>
	EitherT<L, Md, A> EitherT<L, Md, A>::join(const EitherT<L, Md, EitherT<L, Md, A> >& eea)
	{
//		return EitherT<L, Md, A>(eea.run() >> (runEitherT >> Monad<Md>::fmap) >> Monad<Md>::join);
		const auto id = [](const EitherT<L, Md, A>& ea) { return ea; };

		return EitherT<L, Md, A>(eea.run() >> ((either << EitherT<L, Md, A>::makeLeft << id) >> Monad<Md>::fmap) >> (runEitherT >> Monad<Md>::fmap) >> Monad<Md>::join);
	}

}

#endif
