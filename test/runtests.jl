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
    dat, names = cadata("smoke")
    @test dat[1, :] == [4, 2, 3, 2]

    dat, names = cadata("author")
    @test dat[1, 1:3] == [550, 116, 147]
end
