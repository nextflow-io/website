<!-- title=More syntax sugar for Nextflow developer!
date=2020-11-02
type=post
tags=nextflow,dsl
status=published
author=Paolo Di Tommaso
icon=paolo.jpg
~~~~~~ -->

# More syntax sugar for Nextflow developer  

Latest Nextflow version 2020.10.0 is the first stable release running on Groovy 3. 

The first benefit of this chnage is that now Nextflow can be compiled and run on any modern Java virtual machine, 
from Java 8 up to latest Java 15! 

Along with this, the new Groovy runtime brings a bunch of syntax enhancements that can be useful 
in the everyday life of pipeline developers. Let's see them more in detail. 


## Improved not operator 

The `!` (not) operator can now prefix the `in` and `instaceof` keyword. 
This makes a bit more concise writing some conditional expression, for example the following snippet: 

```
if( !(x in list) ) {
  // .. 
} 
else if( !(x instanceof String) ) {
  // .. 
}
```

could be replaced by the following: 

```
list =  [10,20,30]

if( x !in list ) {
   // .. 
}
else if( x !instanceof String ) {
   // .. 
}
```

Again, this is a small syntax change which only makes the code a bit more 
readable. 


## Elvis assignment operator

The elvis assigment operator `?=` allows the assignment of a value only if it was 
not previouvsly assigned (or it evaluets to `null`). Consider the following example:  

```
def opts = [foo: 1]

opts.foo ?= 10 
opts.bar ?= 20

assert opts.foo == 1 
assert opts.bar == 20
```

In this snippet the assignment `opts.foo ?= 10` is ignored because the dictionary `opts` already 
contains a value for the `foo` attribute, while the following is assigned as expected.  

In other words this is a shortcut for the following idiom:  

```
if( some_variable != null ) {
  some_variable = 'Hello'
}
```

### Java style lambda expression

Groovy 3 supports the syntax for Java lambda expression. If you don't know what is a Java lamda expression 
don't worry, it's a concept very similar to a Groovy closure, thought with slight differences
both in the syntax and the semantic (i na few words a Groovy closure can modify a variable in the outside scope,
while a Java lamda cannot). 

In terms of syntax a Groovy closure is defined as: 

```
{ it -> SOME_EXPRESSION_HERE }
```

While Java lamba looks like: 

```
it -> { SOME_EXPRESSION_HERE }
```

which can be simplified to the following form when the expression is a single statement: 

```
it -> SOME_EXPRESSION_HERE
```

The good news is that the two syntax are interoperable in many cases and we can use the *lamda*
syntax to get rid-off of the curly bracket parentheses used by the Groovy notation and make our Nextflow 
script more readable. 

For example the following Nextlow idiom: 

```
Channel
    .of( 1,2,3 ) 
    .map { it * it +1 }
    .view { "the value is $it" }
```

Can be rewritting using the labda syntax as: 

```
Channel
    .of( 1,2,3 ) 
    .map( it -> it * it +1 )
    .view( it -> "the value is $it" )
```

Which is a bit more consistent. Note however that the `it ->` implicit argument is now mandatory (while using the closure syntax it could be omitted) and when the operator argument is not *single* value, the lambda requires the 
round parentheses to define the argument e.g. 

```
Channel
    .of( 1,2,3 ) 
    .map( it -> tuple(it * it,  it+1) )
    .view( (a,b) -> "the values are $a and $b" )
```


### Fully support Java streams API  

Java since version 8 provide a [stream library](https://winterbe.com/posts/2014/07/31/java8-stream-tutorial-examples/) that is very powerful and implements some concepts and operators similar to Nextflow channels. 

The main differences between the two is the Nextflow channels and the corresponding operators are *non-blocking* 
i.e. their evaluation is performed asynchronously without blocking your program execution, while Java streams are 
executed in a synchronus manner (at least by default).

A Java stream looks like the following: 

```
assert (1..10).stream()
                .filter(e -> e % 2 == 0)
                .map(e -> e * 2)
                .toList() == [4, 8, 12, 16, 20]

```

Note, in the above example 
[filter](https://docs.oracle.com/javase/8/docs/api/java/util/stream/Stream.html#filter-java.util.function.Predicate-), 
[map](https://docs.oracle.com/javase/8/docs/api/java/util/stream/Stream.html#map-java.util.function.Function-) and 
[toList](https://docs.oracle.com/javase/8/docs/api/java/util/stream/Collectors.html#toList--) 
methods are Java stream operator not the 
[Nextflow](https://www.nextflow.io/docs/latest/operator.html#filter) 
[homonymous](https://www.nextflow.io/docs/latest/operator.html#map) 
[ones](https://www.nextflow.io/docs/latest/operator.html#tolist).


### Java style method reference  

The new runtime allows also the use of the `::` operator to reference an object method. 
This can be useful to pass a method as an argument to a Nextflow operator in a similar manner 
how it was already possible using a closure. For example:

```
Channel
 .of( 'a', 'b', 'c')
 .view( String::toUpperCase )
 ```

The above prints: 

```
  A
  B
  C
```

Because to [view](https://www.nextflow.io/docs/latest/operator.html#filter) operator applied 
the method [toUpperCase](https://docs.oracle.com/javase/8/docs/api/java/lang/String.html#toUpperCase--) 
to each element emitted by the channel. 


### Conclusion 

The new Groovy runtime brings a lot syntax sugar for Nextflow pipelines and allow the use of modern Java 
runtime which deliver better performances and resources usage. 

The ones listed above are only some of them which may be usefull to everyday Nextflow developers. 
If you are curious to learn more about all the changes in the new Groovy parser you can find a 
[this link](https://groovy-lang.org/releasenotes/groovy-3.0.html). 

Finally, a big thanks to the Groovy community for the big effort of developing and maintaining this 
awesome programming environment. 

