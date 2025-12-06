% This function is used to calculate the trait matching between two species

function v_ol = trait_match(z1,z2,sigma)
    v_ol = exp(-(z1-z2)^2/(2*sigma^2));
end