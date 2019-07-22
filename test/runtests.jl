using Ca, Test, LinearAlgebra, CSV

@testset "foo" begin
    x, y = 5, 7
    @test foo(x, y) == 7
    x = "blah"
    # there should be an error when this is called
    @test_throws MethodError foo(x, y)
end

@testset "bar" begin
    z = 4.
    @test bar(z) == 1.
end

@testset "data" begin
    dat, rnames, cnames = cadata("smoke")
    @test dat[1, :] == [4, 2, 3, 2]

    dat, rnames, cnames = cadata("author")
    @test dat[1, 1:3] == [550, 116, 147]
end

@testset "ca" begin
    @test ca(ones(5,5), 2)[1] ≈ 0.8944271909999154

    smoke, rnames, cnames = cadata("smoke")
    @test ca(smoke, 2)[1] ≈ -0.4621229283237114

    @test_throws LinearAlgebra.LAPACKException ca(zeros(5,5), 1)
end
